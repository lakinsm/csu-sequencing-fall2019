## Utility functions that can be optionally used later for producing
## graphs and performing reshaping operations

require(data.table)
require(vegan)

# Call variables from parent scope
`..` <- function (..., .env = sys.parent(2)) {
    get(deparse(substitute(...)), env = .env)
}

# Misc reshape function for data table
melt_dt <- function(D, level_id) {
    temp <- melt(D, variable.name='Sample', value.name='Normalized_Count')
    names(temp) <- c('Name', 'ID', 'Normalized_Count')
    temp <- data.table(cbind(rep(level_id, nrow(temp)), temp))
    names(temp)[1] <- 'Level_ID'
    return(temp)
}

# Calculate the points for convex hulls around data for ordination plots
meg_find_hulls <- function(x) x[chull(x$Ord1, x$Ord2),]

# Function that returns species count, rarefied species count, and alpha diversity measures
# for each sample in the m x n matrix, m = features, n = samples
alpha_rarefaction <- function(X, minlevel, method='invsimpson') {
    S <- specnumber(X, MARGIN=2)
    raremax <- min(colSums(X))
    if( raremax < minlevel ) raremax <- minlevel
    Srare <- rarefy(X, raremax, MARGIN=2)
    Xrare <- t(rrarefy(t(X), raremax))
    alphadiv <- diversity(Xrare, index=method, MARGIN=2)
    return(list(raw_species_abundance=S,
                rarefied_species_abundance=Srare,
                rarefied_data=Xrare,
                alphadiv=alphadiv))
}


heatmap_select_top_counts <- function(X, group_var, sample_var, n) {
    return(X[, tail(.SD, n), by=c(group_var, sample_var)])
}

bar_select_top_counts <- function(X, group_var, n) {
    return(X[, tail(.SD, n), by=group_var])
}


make_sparse <- function(df, rownames, excludes, filter_min_threshold=0.15) {
    local_df <- df[, .SD, .SDcols=!excludes]
    filter_threshold <- quantile(rowSums(local_df), 0.15)
    if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
    chosen <- which(rowSums(local_df) >= filter_threshold )
    ret <- as.data.frame(local_df[chosen, ])
    rownames(ret) <- df[[rownames]][chosen]
    return(ret)
}

data_subset <- function(data_obj, subsets) {
    local_meta <- data.table(pData(data_obj))
    local_subset <- c()          
    for( c in 1:length(subsets) ) {
        
        conditional_terms <- unlist(strsplit(subsets[[c]], ' '))
        conditional_string <- paste('local_meta[[\'', conditional_terms[1],
                                    '\']] ', conditional_terms[2],
                                    ' \'', conditional_terms[3], '\'',
                                    sep='', collapse='')
        if(length(local_subset) > 0) {
            local_subset <- intersect(local_subset, which(eval(parse(text=conditional_string))))
        }
        else {
            local_subset <- which(eval(parse(text=conditional_string)))
        }
    }
    return(data_obj[, local_subset])
}

data_subset_long <- function(data_obj, subsets) {
    local_subset <- c()          
    for( c in 1:length(subsets) ) {
        
        conditional_terms <- unlist(strsplit(subsets[[c]], ' '))
        conditional_string <- paste('data_obj[[\'', conditional_terms[1],
                                    '\']] ', conditional_terms[2],
                                    ' \'', conditional_terms[3], '\'',
                                    sep='', collapse='')
        if(length(local_subset) > 0) {
            local_subset <- intersect(local_subset, which(eval(parse(text=conditional_string))))
        }
        else {
            local_subset <- which(eval(parse(text=conditional_string)))
        }
    }
    return(data_obj[local_subset, ])
}
