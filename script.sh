names="
vol0000
vol0001
vol0002
vol0003
vol0004
vol0005
vol0006
vol0007
vol0008
vol0009
vol0010
vol0011
vol0012
vol0013
vol0014
vol0015
vol0016
vol0017
vol0018
vol0019
vol0020"

for i in $names; do 
mri_vol2surf --src ${i}.nii.gz --srcreg fsaverage/mri/transforms/reg.mni152.2mm.dat --hemi lh --out lh.${i}.mgh --cortex --noreshape --trgsubject fsaverage 
mri_vol2surf --src ${i}.nii.gz --srcreg fsaverage/mri/transforms/reg.mni152.2mm.dat --hemi rh --out rh.${i}.mgh --cortex --noreshape --trgsubject fsaverage 
done

cd /media/sf_Canada_2020/UK_biobank/data/ICA_d25/Yeo7Mapping

flirt -in Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask -ref /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz -out Yeo7_liberal_2mm.nii.gz -applyisoxfm 2



for i in $names; do
fslmaths $i -thr 3 $i
flirt -in $i -ref /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz -out ${i}_1mm.nii.gz -applyisoxfm 1
fslmaths rYeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask -mas ${i}_1mm.nii.gz tmp
for n in `seq 1 7`; do
m=0.5; t=$(echo "($n-$m)"| bc -l)
ut=$(echo "($n+$m)"| bc -l)
fslmaths tmp -thr $t -uthr $ut tmp2
fslstats tmp2 -V >> ${i}.txt
done
#fslstats tmp -V >> ${i}.txt
done
rm tmp.nii.gz tmp2.nii.gz


for i in `ls vol*.nii.gz`; do fslstats $i -C; done > centroid_coords.txt

