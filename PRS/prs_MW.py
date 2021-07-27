import numpy as np, os, pandas as pd
from io import StringIO
from itertools import product
from subprocess import run, CalledProcessError

def get_cohort_SNPs(plink_files, is_plink2):
    if is_plink2:
        try:
            cohort_SNPs = pd.concat([pd.read_table(
                f'{plink_file}.pvar', usecols=['ID', 'REF', 'ALT'],
                index_col='ID') for plink_file in plink_files]).rename(
                columns=dict(REF='A2', ALT='A1'))
        except ValueError:
            cohort_SNPs = pd.concat([pd.read_table(
                StringIO(run(f'grep -v ^## {plink_file}.pvar',
                             check=True, shell=True, capture_output=True
                             ).stdout.decode()),
                usecols=['ID', 'REF', 'ALT'],
                index_col='ID') for plink_file in plink_files]).rename(
                columns=dict(REF='A2', ALT='A1'))
    else:
        cohort_SNPs = pd.concat([pd.read_table(
            f'{plink_file}.bim', usecols=[1, 4, 5], index_col=0)  # ID, ALT, REF
            for plink_file in plink_files])
        cohort_SNPs.columns = 'A1', 'A2'
    # Remove multi-allelic
    cohort_SNPs = cohort_SNPs[~cohort_SNPs.index.duplicated(keep=False)]
    return cohort_SNPs

def get_match_and_flip_combos():
    # Get which allele combos (sumstats ref, sumstats alt, 1000G ref, 1000G alt)
    # match, and which should be flipped
    # Adapted from https://github.com/bulik/ldsc/blob/master/ldscore/sumstats.py
    # complementary bases
    COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # bases
    BASES = COMPLEMENT.keys()
    # true iff strand ambiguous
    STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
                        for x in product(BASES, BASES)
                        if x[0] != x[1]}
    # SNPS we want to keep (pairs of alleles)
    VALID_SNPS = {x for x in map(lambda y: ''.join(y), product(BASES, BASES))
                  if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}
    # T iff SNP 1 has the same alleles as SNP 2
    # (allowing for strand or ref allele flip).
    MATCH_ALLELES = {x for x in map(lambda y: ''.join(y),
                                    product(VALID_SNPS, VALID_SNPS))
                     # strand and ref match
                     if ((x[0] == x[2]) and (x[1] == x[3])) or
                     # ref match, strand flip
                     ((x[0] == COMPLEMENT[x[2]]) and (
                             x[1] == COMPLEMENT[x[3]])) or
                     # ref flip, strand match
                     ((x[0] == x[3]) and (x[1] == x[2])) or
                     ((x[0] == COMPLEMENT[x[3]]) and
                      (x[1] == COMPLEMENT[x[2]]))}  # strand and ref flip
    # T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
    FLIP_ALLELES = {''.join(x) for x in MATCH_ALLELES
                    # strand match
                    if ((x[0] == x[3]) and (x[1] == x[2])) or
                    # strand flip
                    ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}
    return MATCH_ALLELES, FLIP_ALLELES

MATCH_ALLELES, FLIP_ALLELES = get_match_and_flip_combos()

def get_prs(phenotype_name, sumstats_file, prs_weight_file, prs_score_file,
            plink_files, p_thresh=0.05, is_plink2=True):
    # See https://choishingwan.github.io/PRS-Tutorial,
    # https://www.sciencedirect.com/science/article/pii/S0002929719304227 and
    # https://ars.els-cdn.com/content/image/1-s2.0-S0002929719304227-mmc1.pdf
    if isinstance(plink_files, str):
        plink_files = plink_files,
    else:
        plink_files = tuple(plink_files)
    # 0. Return cached PRSs if they exist, or print what we're generating
    if os.path.exists(prs_score_file):
        PRS = pd.read_table(prs_score_file, index_col=0, header=None,
                            squeeze=True)
        return PRS
    if os.path.exists(prs_weight_file):
        print(f'Generating {prs_score_file} using cached {prs_weight_file}...')
    else:
        print(f'Generating {prs_weight_file} and {prs_score_file}...')
        os.makedirs('PRS', exist_ok=True)
        # 1. Get column indices
        for sep in '\t ,':
            columns = pd.read_csv(sumstats_file, sep=sep, nrows=0).columns
            if len(columns) != 1: break
        columns = columns.str.upper()  # case-insensitive
        SNP_column = np.flatnonzero(columns.isin([
            'SNP', 'ID', 'SNPID', 'MARKERNAME']))[0] + 1
        A1_column = np.flatnonzero(columns.isin(['A1', 'REF']))[0] + 1
        A2_column = np.flatnonzero(columns.isin(['A2', 'ALT']))[0] + 1
        beta_column = np.flatnonzero(columns.isin([
            'OR', 'BETA', 'LOGOR']))[0] + 1
        p_column = np.flatnonzero(columns.isin(['P', 'PVAL']))[0] + 1
        # 2. Filter to p < threshold with awk, and subset to SNP,A1,A2,BETA,P
        print('Filtering...')
        awk_sumstats_file = f'<(zcat {sumstats_file})' \
            if sumstats_file.endswith('.gz') else sumstats_file
        filtered_sumstats_file = f'PRS/{phenotype_name}.filtered_sumstats.txt'
        effect_column_string = f'log(${beta_column})' \
            if 'OR' in columns else f'${beta_column}'
        run(f"awk -F '{sep}' -v OFS='\t' '${p_column} < {p_thresh} "
            f"{{print ${SNP_column},toupper(${A1_column}),"
            f"toupper(${A2_column}),{effect_column_string},${p_column}}}' "
            f"{awk_sumstats_file} > {filtered_sumstats_file}",
            check=True, shell=True, executable='/bin/bash')
        # 3. Harmonize alleles
        print('Harmonizing...')
        sumstats = pd.read_table(filtered_sumstats_file,
                                 names=['SNP', 'A1', 'A2', 'BETA', 'P'])
        assert len(sumstats) > 0
        sumstats.set_index('SNP', inplace=True)
        sumstats = sumstats[~sumstats.index.duplicated(keep=False)]  # multi-all
        cohort_SNPs = get_cohort_SNPs(plink_files, is_plink2)
        SNPs_in_both = sumstats.index.intersection(cohort_SNPs.index)
        assert len(SNPs_in_both) > 0, \
            f'{sumstats.index[0]=}, {cohort_SNPs.index[0]=}'
        sumstats = sumstats.loc[SNPs_in_both]
        cohort_SNPs = cohort_SNPs.loc[SNPs_in_both]
        quadruplets = sumstats.A1 + sumstats.A2 + \
                      cohort_SNPs.A1 + cohort_SNPs.A2
        keep = quadruplets.isin(MATCH_ALLELES)  # also removes indels
        sumstats = sumstats[keep]
        quadruplets = quadruplets[keep]
        assert len(keep) > 0, f'{quadruplets.iloc[0]=}'
        print(f'Kept {100 * keep.mean():.1f}% of SNPs that are non-ambiguous.')
        flip = quadruplets.isin(FLIP_ALLELES)
        sumstats.BETA[flip] = -sumstats.BETA[flip]
        sumstats['A1'] = cohort_SNPs.A1
        print(f'Flipped {100 * flip.mean():.1f}% of SNPs to match cohort.')
        # 4. LD-prune to r^2 < 0.5 with plink
        # (note: this implicitly filters out SNPs that are not in cohort)
        print('LD-pruning...')
        sumstats.index.to_frame().to_csv(
            f'PRS/{phenotype_name}.SNPs.txt', index=False, header=False)
        for file_index, plink_file in enumerate(plink_files):
            run(f'rm -f PRS/{phenotype_name}.{file_index}.prune.in',
                shell=True, check=True)
            try:
                run(f'plink2 '
                    f'--seed 0 '
                    f'{"--pfile" if is_plink2 else "--bfile"} {plink_file} '
                    f'--extract PRS/{phenotype_name}.SNPs.txt '
                    f'--indep-pairwise 500kb 0.5 '
                    f'--out PRS/{phenotype_name}.{file_index}',
                    shell=True, check=True, executable='/bin/bash')
            except CalledProcessError:
                print(f'WARNING: no variants from {plink_file} '
                      f'remain; continuing')
                continue
        # Filter sumstats to these SNPs, and save rs/A1/beta as PRS weight file
        pruned_SNPs = pd.concat([pd.read_table(
            filename, header=None, squeeze=True)
            for file_index in range(len(plink_files))
            if os.path.isfile(
                filename := f'PRS/{phenotype_name}.{file_index}.prune.in')
               and os.path.getsize(filename) > 0])
        sumstats = sumstats[sumstats.index.isin(pruned_SNPs)]
        os.makedirs(os.path.dirname(prs_weight_file), exist_ok=True)
        sumstats[['A1', 'BETA']].to_csv(prs_weight_file, sep='\t', header=False)
    # 5. Score PRS on cohort
    print('Scoring...')
    file_PRSs = {}
    for file_index, plink_file in enumerate(plink_files):
        try:
            run(f'plink2 '
                f'{"--pfile" if is_plink2 else "--bfile"} {plink_file} '
                f'--extract <(cut -f1 {prs_weight_file}) '
                f'--score {prs_weight_file} header variance-standardize '
                f'cols=scoresums '
                f'--out PRS/{phenotype_name}.{file_index}',
                shell=True, check=True, executable='/bin/bash')
        except CalledProcessError:
            print(f'WARNING: no variants from {plink_file} remain; continuing')
            continue
        # Read PRS
        file_PRS = pd.read_table(
            f'PRS/{phenotype_name}.{file_index}.sscore', delim_whitespace=True,
            usecols=['#IID', 'SCORE1_SUM'], index_col='#IID', squeeze=True)
        file_PRSs[plink_file] = file_PRS
    # 6. Sum PRS across files
    print('Aggregating...')
    file_PRSs = pd.DataFrame(file_PRSs, columns=file_PRSs)
    assert not file_PRSs.isna().values.any()
    # noinspection PyArgumentList
    PRS = file_PRSs.sum(axis=1)
    # 7. Save
    os.makedirs(os.path.dirname(prs_score_file), exist_ok=True)
    PRS.to_csv(prs_score_file, sep='\t', header=False)
    return PRS

if __name__ == '__main__':
    # Note: sumstats files must contain the following columns (case doesn't
    # matter): SNP/ID/MARKERNAME (rs number), A1, A2, OR/BETA/LOGOR, P/PVAL
    sumstats_files = {
        'MDD': 'daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz',
        'BIP': 'daner_PGC_BIP32b_mds7a_0416a.gz',
        'SCZ': 'PGC3_SCZ_wave3_public.v2.tsv.gz',
        'ANX': 'daner_woautism_ad_sd8-sd6_woautismad_cleaned.gz',
        'PTSD': 'psych_sumstats/pts_eur_freeze2_overall.results.gz'}
    plink_files = [f'ukb/genomic_data_imp/pgen_white_all/chr{chrom}'
                   for chrom in range(1, 23)]
    # Generate PRSs
    for phenotype_name, sumstats_file in sumstats_files.items():
        p_threshs = {5e-5: 'PRS/5e-5', 5e-4: 'PRS/5e-4', 5e-3: 'PRS/5e-3',
                     0.05: 'PRS/0.05', 0.5: 'PRS/0.5'}
        for p_thresh, prs_dir in p_threshs.items():
            prs_weight_file = f'{prs_dir}/{phenotype_name}.prs.weights'
            prs_score_file = f'{prs_dir}/{phenotype_name}.prs.scores'
            generate_PRS = not os.path.exists(prs_score_file)
            if generate_PRS: print(f'p = {p_thresh}, {phenotype_name}:')
            PRS = get_prs(phenotype_name, sumstats_file, prs_weight_file,
                          prs_score_file, plink_files, p_thresh)

