library(data.table)
library(ggplot2)
library(biomformat)

set.seed(2718)

setwd('/mnt/datasets/coursework/sequencing/')

dat = read_biom('mothur_output.biom')
obs = biom_data(dat)
features = unlist(lapply(dat$rows, function(x){
  return(paste0(x$metadata$taxonomy, collapse='|'))
}))
samples = unlist(lapply(dat$columns, function(x){
  return(x$id)
}))
obs_dt = data.table(as.matrix(obs))
obs_dt[, Feature := features]
obs_dt = obs_dt[, lapply(.SD, sum), by='Feature']
obs_matrix = as.matrix(obs_dt[, .SD, .SDcols=!'Feature'])
rownames(obs_matrix) = obs_dt$Feature

source('scripts/utility.R')
source('scripts/graphing.R')
source('scripts/statistical_analysis.R')

metadata = data.table(
  SampleName=colnames(obs_matrix),
  Group=c('Soil 1', 'Soil 1', 'Soil 1',
          'Soil 2', 'Soil 2', 'Soil 2',
          'Soil 3', 'Soil 4', 'Soil 4',
          'Water', 'Zymo Mock')
)

lib_sizes = data.table(
  SampleName=colnames(obs_matrix),
  RemovedReads=c(102336, 16, 130528, 140880, 124, 140176, 23800, 162136,
                222184, 8008, 102024) / 4
)
lib_sizes[, AnalyticData := colSums(obs_matrix)]
lib_sizes$RemovedReads = lib_sizes$RemovedReads - lib_sizes$AnalyticData

obs_mr = newMRexperiment(obs_matrix[rowSums(obs_matrix) > 0, ])  # remove zero-count features

cumNorm(obs_mr)  # Cumulative Sum Scaling

obs_norm = data.table(MRcounts(obs_mr, norm=T))  # Extract normalized data
obs_raw = data.table(MRcounts(obs_mr, norm=F))  # Save raw data too

# Remove features below the 5th percentile of counts
obs_norm_filt = obs_norm[rowSums(obs_norm) >= quantile(rowSums(obs_norm), 0.05), ]
obs_raw_filt = obs_raw[rowSums(obs_raw) >= quantile(rowSums(obs_raw), 0.05), ]

# Add back in the feature names into the id column (minus those we removed in the previous step)
obs_norm_filt[, ID := rownames(obs_matrix)[rowSums(obs_norm) >= quantile(rowSums(obs_norm), 0.05)]]
obs_raw_filt[, ID := rownames(obs_matrix)[rowSums(obs_raw) >= quantile(rowSums(obs_raw), 0.05)]]

# Save the feature names in a separate variable for later
obs_taxonomy = data.table(ID=rownames(obs_matrix))

# Create columns for all taxonomic ranks
obs_taxonomy = setDT(obs_taxonomy)[, c('Domain',
                        'Phylum',
                        'Class',
                        'Order',
                        'Family',
                        'Genus') := tstrsplit(ID, '|', type.convert=T, fixed=T)]

# Set key similar to a unique key in a SQL table.  Helps with table joins.
setkey(obs_taxonomy, ID)
setkey(obs_norm_filt, ID)
setkey(obs_raw_filt, ID)

# Left SQL join (keep only features in the count table, but pull in the taxonomic columns from taxonomy)
obs_norm_filt = obs_taxonomy[obs_norm_filt]
obs_raw_filt = obs_taxonomy[obs_raw_filt]


# Split data into taxonomic levels, keeping only relevant information
# Do this for both raw and normalized data tables
# obs_norm_filt_domain = obs_norm_filt[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:7]
# obs_norm_filt_domain_analytic = newMRexperiment(counts=obs_norm_filt_domain[, .SD, .SDcols=!'Domain'])
# rownames(obs_norm_filt_domain_analytic) = obs_norm_filt_domain$Domain

obs_norm_filt_phylum = obs_norm_filt[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:7]
obs_norm_filt_phylum_analytic = newMRexperiment(counts=obs_norm_filt_phylum[, .SD, .SDcols=!'Phylum'])
rownames(obs_norm_filt_phylum_analytic) = obs_norm_filt_phylum$Phylum

obs_norm_filt_class = obs_norm_filt[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:7]
obs_norm_filt_class_analytic = newMRexperiment(counts=obs_norm_filt_class[, .SD, .SDcols=!'Class'])
rownames(obs_norm_filt_class_analytic) = obs_norm_filt_class$Class

obs_norm_filt_order = obs_norm_filt[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:7]
obs_norm_filt_order_analytic = newMRexperiment(counts=obs_norm_filt_order[, .SD, .SDcols=!'Order'])
rownames(obs_norm_filt_order_analytic) = obs_norm_filt_order$Order

obs_norm_filt_family = obs_norm_filt[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:7]
obs_norm_filt_family_analytic = newMRexperiment(counts=obs_norm_filt_family[, .SD, .SDcols=!'Family'])
rownames(obs_norm_filt_family_analytic) = obs_norm_filt_family$Family

obs_norm_filt_genus = obs_norm_filt[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:7]
obs_norm_filt_genus_analytic = newMRexperiment(counts=obs_norm_filt_genus[, .SD, .SDcols=!'Genus'])
rownames(obs_norm_filt_genus_analytic) = obs_norm_filt_genus$Genus
# 
# obs_raw_filt_domain = obs_raw_filt[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:7]
# obs_raw_filt_domain_analytic = newMRexperiment(counts=obs_raw_filt_domain[, .SD, .SDcols=!'Domain'])
# rownames(obs_raw_filt_domain_analytic) = obs_raw_filt_domain$Domain

obs_raw_filt_phylum = obs_raw_filt[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:7]
obs_raw_filt_phylum_analytic = newMRexperiment(counts=obs_raw_filt_phylum[, .SD, .SDcols=!'Phylum'])
rownames(obs_raw_filt_phylum_analytic) = obs_raw_filt_phylum$Phylum

obs_raw_filt_class = obs_raw_filt[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:7]
obs_raw_filt_class_analytic = newMRexperiment(counts=obs_raw_filt_class[, .SD, .SDcols=!'Class'])
rownames(obs_raw_filt_class_analytic) = obs_raw_filt_class$Class

obs_raw_filt_order = obs_raw_filt[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:7]
obs_raw_filt_order_analytic = newMRexperiment(counts=obs_raw_filt_order[, .SD, .SDcols=!'Order'])
rownames(obs_raw_filt_order_analytic) = obs_raw_filt_order$Order

obs_raw_filt_family = obs_raw_filt[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:7]
obs_raw_filt_family_analytic = newMRexperiment(counts=obs_raw_filt_family[, .SD, .SDcols=!'Family'])
rownames(obs_raw_filt_family_analytic) = obs_raw_filt_family$Family

obs_raw_filt_genus = obs_raw_filt[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:7]
obs_raw_filt_genus_analytic = newMRexperiment(counts=obs_raw_filt_genus[, .SD, .SDcols=!'Genus'])
rownames(obs_raw_filt_genus_analytic) = obs_raw_filt_genus$Genus

# Combine data tables into a list data structure for easier looping/manipulation
# We melt (wide-to-long reformat) these data tables to work with heatmaps and barcharts.
obs_norm_melted_analytic = rbind(melt_dt(MRcounts(obs_norm_filt_phylum_analytic), 'Phylum'),
                                 melt_dt(MRcounts(obs_norm_filt_class_analytic), 'Class'),
                                 melt_dt(MRcounts(obs_norm_filt_order_analytic), 'Order'),
                                 melt_dt(MRcounts(obs_norm_filt_family_analytic), 'Family'),
                                 melt_dt(MRcounts(obs_norm_filt_genus_analytic), 'Genus'))
obs_raw_melted_analytic = rbind(melt_dt(MRcounts(obs_raw_filt_phylum_analytic), 'Phylum'),
                                melt_dt(MRcounts(obs_raw_filt_class_analytic), 'Class'),
                                melt_dt(MRcounts(obs_raw_filt_order_analytic), 'Order'),
                                melt_dt(MRcounts(obs_raw_filt_family_analytic), 'Family'),
                                melt_dt(MRcounts(obs_raw_filt_genus_analytic), 'Genus'))

# These are the wide (unmelted) data tables
obs_norm_analytic = c(obs_norm_filt_phylum_analytic,
                      obs_norm_filt_class_analytic,
                      obs_norm_filt_order_analytic,
                      obs_norm_filt_family_analytic,
                      obs_norm_filt_genus_analytic)
obs_raw_analytic = c(obs_raw_filt_phylum_analytic,
                     obs_raw_filt_class_analytic,
                     obs_raw_filt_order_analytic,
                     obs_raw_filt_family_analytic,
                     obs_raw_filt_genus_analytic)


# Name the data lists for easier access
names(obs_norm_analytic) = c('Phylum', 'Class', 'Order', 'Family', 'Genus')
names(obs_raw_analytic) = c('Phylum', 'Class', 'Order', 'Family', 'Genus')

# Make sure the column names are the same order for both the metadata and data frames.
# Add back in the feature names where needed
for(l in 1:length(obs_norm_analytic)) {
  sample_idx = match(colnames(MRcounts(obs_norm_analytic[[l]])), metadata$SampleName)
  pData(obs_norm_analytic[[l]]) = data.frame(metadata[sample_idx, .SD, .SDcols=!'SampleName'])
  rownames(pData(obs_norm_analytic[[l]])) = metadata[['SampleName']][sample_idx]
  fData(obs_norm_analytic[[l]]) = data.frame(Feature=rownames(MRcounts(obs_norm_analytic[[l]])))
  rownames(fData(obs_norm_analytic[[l]])) = rownames(MRcounts(obs_norm_analytic[[l]]))
}

for(l in 1:length(obs_raw_analytic)) {
  sample_idx = match(colnames(MRcounts(obs_raw_analytic[[l]])), metadata$SampleName)
  pData(obs_raw_analytic[[l]]) = data.frame(metadata[sample_idx, .SD, .SDcols=!'SampleName'])
  rownames(pData(obs_raw_analytic[[l]])) = metadata[['SampleName']][sample_idx]
  fData(obs_raw_analytic[[l]]) = data.frame(Feature=rownames(MRcounts(obs_raw_analytic[[l]])))
  rownames(fData(obs_raw_analytic[[l]])) = rownames(MRcounts(obs_raw_analytic[[l]]))
}

# Same thing for the metadata table, then set the SQL key as SampleName
metadata <- data.table(metadata[match(colnames(obs_matrix), metadata[['SampleName']])])
setkeyv(metadata, 'SampleName')

# The following functions are custom and are defined in the other R files in the Scripts folder.
## Graphing
# Alpha rarefaction
meg_alpha_rarefaction(obs_raw_analytic, # Data table list (wide format)
                      names(obs_raw_analytic), # The list of taxonomic level names
                      metadata, # Metadata table
                      'SampleName', # Column name where sample ID's are stored
                      'Group', # Factor to analyze (color)
                      list(), # Any subsetting rules to apply (usually left as an empty list)
                      'outputs/raw/', # Output folder
                      'Microbiome') # Analysis to run (choice of Microbiome or Resistome)


# PCA (arguments are the same as above, with changes noted below)
meg_ordination(obs_norm_analytic,
               names(obs_norm_analytic),
               metadata,
               'SampleName',
               'Group',
               list(),
               'outputs/normalized/',
               'Microbiome',
               'PCA')  # The analysis type (choice of PCA or NMDS)

meg_ordination(obs_raw_analytic,
               names(obs_raw_analytic),
               metadata,
               'SampleName',
               'Group',
               list(),
               'outputs/raw/',
               'Microbiome',
               'PCA')


# NMDS
meg_ordination(obs_norm_analytic,
               names(obs_norm_analytic),
               metadata,
               'SampleName',
               'Group',
               list(),
               'outputs/normalized/',
               'Microbiome',
               'NMDS')

meg_ordination(obs_raw_analytic,
               names(obs_raw_analytic),
               metadata,
               'SampleName',
               'Group',
               list(),
               'outputs/raw/',
               'Microbiome',
               'NMDS')


# Heatmaps (Note we are looping over the data lists here, unlike the previous functions)
for(l in 1:length(obs_norm_analytic)) {
  meg_heatmap(obs_norm_melted_analytic,  # Note that this is the melted data table list
              metadata,
              'SampleName',
              'Group',
              names(obs_norm_analytic)[l],
              list(),
              'outputs/normalized/',
              'Microbiome'
  )
}

for(l in 1:length(obs_raw_analytic)) {
  meg_heatmap(obs_raw_melted_analytic,
              metadata,
              'SampleName',
              'Group',
              names(obs_raw_analytic)[l],
              list(),
              'outputs/raw/',
              'Microbiome'
  )
}


# Barplots (Note we are looping over the data lists here, unlike the previous functions)
for(l in 1:length(obs_norm_analytic)) {
  meg_barplot(obs_norm_melted_analytic,
              metadata,
              'SampleName',
              'Group',
              names(obs_norm_analytic)[l],
              list(),
              'outputs/normalized/',
              'Microbiome')
}

for(l in 1:length(obs_raw_analytic)) {
  meg_barplot(obs_raw_melted_analytic,
              metadata,
              'SampleName',
              'Group',
              names(obs_raw_analytic)[l],
              list(),
              'outputs/raw/',
              'Microbiome')
}

# Output the matrices to a separate folder for each taxonomic level
for(l in 1:length(obs_norm_analytic)) {
  write.csv(MRcounts(obs_norm_analytic[[l]]),
            paste('matrices/normalized/', names(obs_norm_analytic)[l], '.csv', sep='', collapse = ''),
            row.names=T)
}
for(l in 1:length(obs_raw_analytic)) {
  write.csv(MRcounts(obs_raw_analytic[[l]]),
            paste('matrices/raw/', names(obs_raw_analytic)[l], '.csv', sep='', collapse = ''),
            row.names=T)
}


# Alpha rarefaction
library_sizes = colSums(obs_raw_analytic$Genus)
sample_vec = round(seq(10, max(library_sizes), length.out=20))
local_data = MRcounts(obs_raw_analytic$Genus)
rarefaction = data.table(
  SampleName=character(),
  SampleSize=numeric(),
  MeanSpeciesRichness=numeric(),
  SESpeciesRichness=numeric()
)
for(s in sample_vec) {
  bootstrap_matrix = data.table()
  res = rarefy(local_data, sample=s, se=T, MARGIN=2)
  rarefaction = rbind(rarefaction,
                      data.table(
                        SampleName=colnames(res),
                        SampleSize=rep(s, ncol(res)),
                        MeanSpeciesRichness=res['S', ],
                        SESpeciesRichness=res['se', ]
                      ))
}
rarefaction$MeanSpeciesRichness[rarefaction$SESpeciesRichness == 0] = NA

g = ggplot(rarefaction, aes(x=SampleSize, y=MeanSpeciesRichness, 
                            color=SampleName)) +
  geom_line() +
  geom_errorbar(aes(ymin=MeanSpeciesRichness-SESpeciesRichness,
                    ymax=MeanSpeciesRichness+SESpeciesRichness),
                alpha=0.3) +
  theme(
    legend.title=element_text(size=18),
    legend.text=element_text(size=14),
    axis.title=element_text(size=18),
    axis.text.x=element_text(size=14, hjust=1, angle=45),
    axis.text.y=element_text(size=14)
  )
print(g)


# Library size plot
mlib = melt(lib_sizes)
g = ggplot(mlib, aes(x=SampleName, y=value, fill=variable)) +
  geom_bar(stat='identity') +
  geom_text(aes(x=SampleName, y=totals + 2000,
                label=percentiles,
                fill=NULL), data=data.table(
                  SampleName=colnames(obs_matrix),
                  percentiles=paste(round(100*lib_sizes$RemovedReads / (lib_sizes$RemovedReads + lib_sizes$AnalyticData), 1), '%', sep=''),
                  totals=lib_sizes$RemovedReads + lib_sizes$AnalyticData
                ), size=4) +
  theme(
    legend.title=element_text(size=18),
    legend.text=element_text(size=14),
    axis.title=element_text(size=18),
    axis.text.x=element_text(size=14, hjust=1, angle=45),
    axis.text.y=element_text(size=14)
  ) +
  scale_fill_discrete(name='Read Type') +
  ylab('Count of Paired-End Reads') +
  xlab('Sample Name')
print(g)



# Features at each level
feat_count = data.table(
  TaxonomicLevel=names(obs_norm_analytic),
  ClassifiedCounts=unlist(lapply(obs_norm_analytic, function(X){
    selected_features = rownames(X)[!grepl('(unclassified|uncultured|unknown)', rownames(X), ignore.case=T, perl=T)]
    local_matrix = MRcounts(X)
    return(sum(as.vector(local_matrix[selected_features, ])))
  }))
)
feat_count[, PercentTotal := paste(round(100*ClassifiedCounts/max(ClassifiedCounts), 1), '%', sep='') ]


feat_count$TaxonomicLevel = factor(feat_count$TaxonomicLevel,
                                   levels=c('Phylum', 'Class', 'Order', 'Family',
                                            'Genus'), ordered=T)
g = ggplot(feat_count, aes(x=TaxonomicLevel, y=ClassifiedCounts)) +
  geom_bar(stat='identity') +
  geom_text(aes(y=ClassifiedCounts + 5000, label=PercentTotal), size=4) +
  theme(
    legend.title=element_text(size=18),
    legend.text=element_text(size=14),
    axis.title=element_text(size=18),
    axis.text.x=element_text(size=14, hjust=1, angle=45),
    axis.text.y=element_text(size=14)
  ) +
  xlab('Taxonomic Level') +
  ylab('Classified Counts')
print(g)


# Looking at Zymo mock community
zym = MRcounts(obs_norm_analytic$Genus)
zym = as.data.table(zym[, 'Zymo_mock'])
zym[, TaxonomicID := rownames(obs_norm_analytic$Genus)]
zym = zym[V1 > 0, ]
zym = zym[!grepl('unclassified', TaxonomicID, ignore.case=T)]
zym[, NormalizedCount := V1]
zym$V1 = paste(round(100*zym$V1 / sum(zym$V1), 1), '%', sep='')
colnames(zym) = c('Percent', 'TaxonomicID', 'NormalizedCount')

g = ggplot(zym, aes(x=TaxonomicID, y=NormalizedCount, fill=TaxonomicID)) +
  geom_bar(stat='identity') +
  geom_text(aes(y=NormalizedCount + 100, label=Percent), size=4) +
  theme(
    legend.title=element_text(size=18),
    legend.text=element_text(size=14),
    axis.title=element_text(size=18),
    axis.text.x=element_text(size=14, hjust=1, angle=45),
    axis.text.y=element_text(size=14),
    legend.position='None',
    plot.title=element_text(size=20, hjust=0.5)
  ) +
  xlab('Genus') +
  ylab('Normalized Count') +
  ggtitle('Zymo Mock Community')
print(g)


# Water negative control
wat = MRcounts(obs_norm_analytic$Order)
wat = as.data.table(wat[, 'Water'])
wat[, TaxonomicID := rownames(obs_norm_analytic$Order)]
wat = wat[V1 > 0, ]
wat = wat[!grepl('unclassified', TaxonomicID, ignore.case=T)]
wat[, NormalizedCount := V1]
wat$V1 = paste(round(100*wat$V1 / sum(wat$V1), 1), '%', sep='')
colnames(wat) = c('Percent', 'TaxonomicID', 'NormalizedCount')
wat = wat[NormalizedCount > quantile(NormalizedCount, 0.5)]

g = ggplot(wat, aes(x=TaxonomicID, y=NormalizedCount, fill=TaxonomicID)) +
  geom_bar(stat='identity') +
  geom_text(aes(y=NormalizedCount + 100, label=Percent), size=4) +
  theme(
    legend.title=element_text(size=18),
    legend.text=element_text(size=14),
    axis.title=element_text(size=18),
    axis.text.x=element_text(size=14, hjust=1, angle=45),
    axis.text.y=element_text(size=14),
    legend.position='None',
    plot.title=element_text(size=20, hjust=0.5)
  ) +
  xlab('Genus') +
  ylab('Normalized Count') +
  ggtitle('Water Negative Control')
print(g)



