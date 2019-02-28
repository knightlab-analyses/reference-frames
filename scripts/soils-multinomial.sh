##Natbiotech revisions; run Central Park Soil data through Songbird to idenitfy missed/false associations
results_dir=songbird-microbial_biomass+ph+tot_nitro+tot_org_carb+water_content_soil-epoch10000-lr-3-batch10_minf24_deblur
songbird multinomial \
	--formula 'microbial_biomass+ph+tot_nitro+tot_org_carb+water_content_soil' \
	--input-biom ../data/soils/ref_table.biom \
	--metadata-file ../data/soils/metadata.txt \
	--batch-size 20 \
	--epoch 10000 \
    	--learning-rate 1e-3 \
    	--min-feature-count 24 \
	--summary-dir $results_dif

mv $results_dir/beta.csv ../results/soil_differentials.csv
