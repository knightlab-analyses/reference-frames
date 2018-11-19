

songbird multinomial \
  --input-biom ../data/oral_deblurred.biom \
  --metadata-file ../data/oral_metadata.txt \
  --formula 'C(brushing_event)' \
  --training-column Test \
  --epoch 40000 \
  --learning-rate 1e-3 \
  --batch-size 3 \
  --beta-prior 10 \
  --summary-dir ../results/oral-resulrs
