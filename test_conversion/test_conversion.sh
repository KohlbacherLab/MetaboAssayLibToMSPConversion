#!/bin/bash

echo 'test - .pqp'
python ../parseMetaboAssayLibToMsp.py \
-in ./assaylibrary_t3_decoy_test.pqp \
-out ./msp_pqp_test.msp \

echo '-----------------------------'

echo 'test - .TraML'
python ../parseMetaboAssayLibToMsp.py \
-in ./assaylibrary_t3_decoy_test.TraML \
-out ./msp_from_TraML_test.msp \

echo '-----------------------------'

echo 'test - .tsv'
python ../parseMetaboAssayLibToMsp.py \
-in ./assaylibrary_t3_decoy_test.tsv \
-out ./msp_from_tsv_test.msp \

