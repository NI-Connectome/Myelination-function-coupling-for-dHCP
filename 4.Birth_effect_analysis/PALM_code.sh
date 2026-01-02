#!/bin/zs

echo "convert nii to gii"

palm -i /data/MFC_left.shape.gii \
     -s /data/hemi-left_midthickness.5k.surf.gii \
     -i /data/MFC_right.shape.gii \
     -s /data/hemi-right_midthickness.5k.surf.gii \
     -d /data/designmatrix_test.mat \
     -t /data/contrasts_test.con \
     -o /result/MFC \
     -n 5000 \
     -T \
     -tfce2d \
     -precision double \
     -logp \
     -approx gamma \
     -ise
