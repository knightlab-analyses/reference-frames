

songbird multinomial \
  --input-biom ../data/oral_trimmed_deblur.biom \
  --metadata-file ../data/oral_trimmed_metadata.txt \
  --formula 'C(brushing_event)' \
  --training-column Test \
  --epochs 40000 \
  --learning-rate 1e-3 \
  --batch-size 3 \
  --differential-prior 10 \
  --summary-dir ../results/oral-results
