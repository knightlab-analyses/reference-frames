songbird multinomial \
  --input-biom ../data/byrd_skin_table.biom \
  --metadata-file ../data/byrd_metadata.txt \
  --formula "C(Timepoint, Treatment('F'))" \
  --training-column Test \
  --epochs 5000 \
  --learning-rate 1e-3 \
  --batch-size 5 \
  --beta-prior 10 \
  --summary-dir ../results/byrd-results
