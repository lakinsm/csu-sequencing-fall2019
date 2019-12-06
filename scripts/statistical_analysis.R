require(data.table)
require(metagenomeSeq)

meg_fitZig <- function(data_list,
                       data_names,
                       metadata,
                       zero_mod,
                       data_mod,
                       filter_min_threshold,
                       contrast_list,
                       random_effect_var,
                       outdir,
                       analysis_name,
                       analysis_subset,
                       data_type,
                       pval,
                       top_hits) {
    settings <- zigControl(maxit=50, verbose=F)
    
    local_obj <- data_list
    res <- list()
    for( l in 1:length(local_obj) ) {
        
        filter_threshold <- 0
        if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
        local_obj[[l]] <- local_obj[[l]][which(rowSums(MRcounts(local_obj[[l]])) >= filter_threshold ), ]
        
        if(length(analysis_subset) > 0) {
            local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
        }
        
        
        col_selection <- as.integer(which(colSums(MRcounts(local_obj[[l]]) > 0) > 1))
        local_obj[[l]] <- local_obj[[l]][, col_selection]
        
        tryCatch(
            {
                mod_select <- model.matrix(eval(parse(text=data_mod)), data=pData(local_obj[[l]]))
                zero_mod_select <- zero_mod[col_selection, ]
            },
            error=function(e) {
                print(paste('Not enough non-zero features for one contrast group for ', data_type, ' ', data_names[l], ' ', analysis_name,
                            sep='', collapse=''))
            }
        )
        
        
        cumNorm(local_obj[[l]])  # This is a placeholder for metagenomeSeq; we don't actually use these values
        
        tryCatch(
            {
                if( is.na(random_effect_var) ) {
                    res[[l]] <- fitZig(obj=local_obj[[l]],
                                       mod=mod_select,
                                       zeroMod=zero_mod_select,
                                       control=settings,
                                       useCSSoffset=F)
                }
                else {
                    res[[l]] <- fitZig(obj=local_obj[[l]],
                                       mod=mod_select,
                                       zeroMod=zero_mod_select,
                                       control=settings,
                                       useCSSoffset=F,
                                       useMixedModel=T,
                                       block=pData(local_obj[[l]])[, random_effect_var])
                }
            },
            error=function(e) {
                print(paste('Model failed to converge for ', data_type, ' ', data_names[l], ' ', analysis_name,
                            sep='', collapse=''))
            },
            finally={
                if( length(res) != l ) {
                    next
                }
            }
        )
        
        local_contrasts <- contrast_list
        local_contrasts[[length(local_contrasts)+1]] <- res[[l]]$fit$design
        names(local_contrasts)[length(local_contrasts)] <- 'levels'
        
        contrast_matrix <- do.call(makeContrasts, local_contrasts)
        colnames(contrast_matrix) <- make.names(contrast_list)
        
        contrast_fit <- contrasts.fit(res[[l]]$fit, contrast_matrix)        
        contrast_fit <- eBayes(contrast_fit)
        
        stats_results <- data.table(
            Node.Name=character(),
            Contrast=character(),
            logFC=numeric(),
            CI.L=numeric(),
            CI.R=numeric(),
            AveExpr=numeric(),
            t=numeric(),
            P.Value=numeric(),
            adj.P.Val=numeric(),
            B=numeric()
        )
        
        for( c in 1:ncol(contrast_fit$contrasts) ) {
            tophits <- topTable(contrast_fit, p.value=pval, confint=T,
                                number=top_hits, sort.by='AveExpr', coef=c)
            
            if( nrow(tophits) > 0) {
                temp_res <- data.table(
                    Node.Name=rownames(tophits),
                    Contrast=rep(colnames(contrast_fit$contrasts)[c], nrow(tophits))
                )
                temp_res <- cbind(temp_res, tophits)
                stats_results <- rbind(stats_results, temp_res)
            }
            else {
                print(paste('No significant results for', data_type,
                            data_names[l], analysis_name,
                            colnames(contrast_fit$contrasts)[c],
                            sep=' ', collapse=''))
            }
        }
        
        if( nrow(stats_results) > 0 ) {
            write.csv(stats_results,
                      file=paste(outdir, '/', analysis_name, '_', data_type, '_',
                                 data_names[l], '_',
                                 contrast_list[1], '_Model_Contrasts.csv',
                                 sep='', collapse=''),
                      quote=F, row.names=F)
        }
    }
}


zero_inflated_gaussian_regression_workflow <- function(experiment_data_lists,
                                                       original_MRexps,
                                                       sample_column_id,
                                                       statistical_analyses,
                                                       stats_output_dir,
                                                       sig_level,
                                                       top_hits,
                                                       lowpass_filter_threshold) {
    if(!is.na(experiment_data_lists[[1]])) {
        # AMR ZIG Regression
        for( a in 1:length(statistical_analyses) ) {
            meg_fitZig(data_list=experiment_data_lists[[1]],
                       data_names=experiment_data_lists[[5]],
                       metadata=experiment_data_lists[[11]],
                       zero_mod=model.matrix(~1 + log(libSize(original_MRexps[[1]]))),
                       data_mod=statistical_analyses[[a]]$model_matrix,
                       filter_min_threshold=lowpass_filter_threshold,
                       contrast_list=statistical_analyses[[a]]$contrasts,
                       random_effect_var=statistical_analyses[[a]]$random_effect,
                       outdir=paste(stats_output_dir, 'AMR', statistical_analyses[[a]]$name,
                                    sep='/', collapse=''),
                       analysis_name=statistical_analyses[[a]]$name,
                       analysis_subset=statistical_analyses[[a]]$subsets,
                       data_type='AMR',
                       pval=sig_level,
                       top_hits=top_hits)
        }
    }
    
    if(!is.na(experiment_data_lists[[6]])) {
        # Microbiome ZIG Regression
        for( a in 1:length(statistical_analyses) ) {
            meg_fitZig(data_list=experiment_data_lists[[6]],
                       data_names=experiment_data_lists[[10]],
                       metadata=experiment_data_lists[[11]],
                       zero_mod=model.matrix(~1 + log(libSize(original_MRexps[[2]]))),
                       data_mod=statistical_analyses[[a]]$model_matrix,
                       filter_min_threshold=lowpass_filter_threshold,
                       contrast_list=statistical_analyses[[a]]$contrasts,
                       random_effect_var=statistical_analyses[[a]]$random_effect,
                       outdir=paste(stats_output_dir, 'Microbiome', statistical_analyses[[a]]$name,
                                    sep='/', collapse=''),
                       analysis_name=statistical_analyses[[a]]$name,
                       analysis_subset=statistical_analyses[[a]]$subsets,
                       data_type='Microbiome',
                       pval=sig_level,
                       top_hits=top_hits)
        }
    }
}


elastic_net_regression_workflow <- function(experiment_data_lists,
                                            sample_column_id,
                                            statistical_analyses,
                                            stats_output_dir) {
    
}


generate_contrasts <- function(analytic_MRexp,
                               feature) {
    values <- data.table(pData(analytic_MRexp))
    values <- as.character(unique(values[[feature]]))
    combs <- utils::combn(values, 2)
    strs <- apply(combs, 2, function(x) {
        paste(paste(feature, x[1], sep='', collapse=''),
              ' - ',
              paste(feature, x[2], sep='', collapse=''),
              sep='', collapse='')
    })
    return(strs)
}


exploratory_zero_inflated_gaussian_regression <- function(analytic_MRexp,
                                                          metadata,
                                                          annotation_level,
                                                          metadata_feature,
                                                          zero_mod,
                                                          data_mod,
                                                          random_effect_var,
                                                          data_type,
                                                          pval,
                                                          top_hits,
                                                          sort_by) {
    settings <- zigControl(maxit=50, verbose=F)
    
    local_MRexp <- analytic_MRexp
    
    col_selection <- as.integer(which(colSums(MRcounts(local_MRexp) > 0) > 1))
    local_MRexp <- local_MRexp[, col_selection]
    pData
    
    mod_select <- model.matrix(eval(parse(text=data_mod)), data=pData(local_MRexp))
    zero_mod_select <- zero_mod[col_selection, ]
    
    cumNorm(local_MRexp)  # This is a placeholder for metagenomeSeq; we don't actually use these values
    
    res <- NA
    
    tryCatch(
        {
            if( is.na(random_effect_var) ) {
                res <- fitZig(obj=local_MRexp,
                                   mod=mod_select,
                                   zeroMod=zero_mod_select,
                                   control=settings,
                                   useCSSoffset=F)
            }
            else {
                res <- fitZig(obj=local_MRexp,
                                   mod=mod_select,
                                   zeroMod=zero_mod_select,
                                   control=settings,
                                   useCSSoffset=F,
                                   useMixedModel=T,
                                   block=pData(local_MRexp)[, random_effect_var])
            }
        },
        error=function(e) {
            print(paste('Model failed to converge for ', data_type, ' ', annotation_level,
                        sep='', collapse=''))
        }
    )
    
    if(is.na(res)) {
        return(data.table(
            Node.Name=character(),
            Contrast=character(),
            logFC=numeric(),
            CI.L=numeric(),
            CI.R=numeric(),
            AveExpr=numeric(),
            t=numeric(),
            P.Value=numeric(),
            adj.P.Val=numeric(),
            B=numeric()))
    }
    
    contrast_list <- generate_contrasts(local_MRexp,
                                        metadata_feature)
    
    local_contrasts <- as.list(contrast_list)
    local_contrasts[[length(local_contrasts)+1]] <- res$fit$design
    names(local_contrasts)[length(local_contrasts)] <- 'levels'
    
    contrast_matrix <- do.call(makeContrasts, local_contrasts)
    colnames(contrast_matrix) <- make.names(contrast_list)
    
    contrast_fit <- contrasts.fit(res$fit, contrast_matrix)        
    contrast_fit <- eBayes(contrast_fit)
    
    stats_results <- data.table(
        Node.Name=character(),
        Contrast=character(),
        logFC=numeric(),
        CI.L=numeric(),
        CI.R=numeric(),
        AveExpr=numeric(),
        t=numeric(),
        P.Value=numeric(),
        adj.P.Val=numeric(),
        B=numeric()
    )
    
    local_sort <- NA
    # 'P-value', 'Effect Size', 'Abundance', 'F-statistic'
    if(sort_by == 'P-value') {
        local_sort <- 'P'
    }
    else if(sort_by == 'Effect Size') {
        local_sort <- 'M'
    }
    else if(sort_by == 'Abundance') {
        local_sort <- 'A'
    }
    else if(sort_by == 'T-statistic') {
        local_sort <- 't'
    }
    
    for( c in 1:ncol(contrast_fit$contrasts) ) {
        tophits <- topTable(contrast_fit, p.value=pval, confint=T,
                            number=top_hits, sort.by=local_sort, coef=c)
        
        if( nrow(tophits) > 0) {
            temp_res <- data.table(
                Node.Name=rownames(tophits),
                Contrast=rep(colnames(contrast_fit$contrasts)[c], nrow(tophits))
            )
            temp_res <- cbind(temp_res, tophits)
            stats_results <- rbind(stats_results, temp_res)
        }
        else {
            print(paste('No significant results for', data_type,
                        annotation_level, colnames(contrast_fit$contrasts)[c],
                        sep=' ', collapse=''))
        }
    }
    
    return(stats_results)
}


generate_permanova_table <- function(data,
                                     data_type,
                                     metadata,
                                     annotation_level,
                                     metadata_feature,
                                     low_pass_filter_threshold,
                                     norm_method,
                                     sample_depth,
                                     strata,
                                     nested_features,
                                     fixed_features) {
    out_dt <- data.table()
    
    # if(data_type == 'Resistome') {
    #     
    #     # Normalization
    #     if(norm_method == 'Rarefaction') {
    #         amr_data_list <- rarefy_normalize_and_extract(sampling_depth=sample_depth,
    #                                                       filtering_quantile=low_pass_filter_threshold,
    #                                                       amr_MRexp=data[[1]])
    #         amr_data_list[[7]] <- data[[2]]
    #     }
    #     else if(norm_method == 'Cumulative Sum Scaling') {
    #         amr_data_list <- CSS_normalize_and_extract(filtering_quantile=low_pass_filter_threshold,
    #                                                    amr_MRexp=data[[1]])
    #         amr_data_list[[7]] <- data[[2]]
    #     }
    #     
    #     # Aggregation and data set creation
    #     analytic_list <- aggregate_and_filter(amr_data_list, data.table(metadata))
    #     
    #     # Distance matrix generation
    #     distances <- list('euclidean', 'bray', 'jaccard', 'manhatten')
    #     dist_mat <- lapply(distances, function(x) {
    #         vegdist(t(MRcounts(analytic_list[[1]][[amr_lookup[annotation_level]]])),
    #                 method = x)
    #     })
    #     
    #     temp <<- dist_mat
    #     # PERMANOVA table generation
    #     strata_string <- ''
    #     if(strata != 'None') {
    #         if(!is.na(nested_features)) {
    #             if(metadata_feature %in% nested_features) {
    #                 strata_string <- paste(list(strata), nested_features, sep=' / ', collapse=' + ')
    #             }
    #         }
    #     }
    #     if(strata_string != '') {
    #         if(length(fixed_features) > 0) {
    #             model_string <- paste(strata_string)
    #         }
    #         else {
    #             model_string <- strata_string
    #         }
    #     }
    #     else {
    #         model_string <- paste(fixed_features, collapse=' + ')
    #     }
    #     res <- lapply(dist_mat, function(x) {
    #         local_model <- paste(' ~', model_string, sep=' ')
    #         if(strata != 'None') {
    #             local_res <- adonis()
    #         }
    #         else {
    #             local_res <- adonis()
    #         }
    #         
    #     })
    #     
    # }
    # else if(data_type == 'Microbiome') {
    #     
    #     # Normalization
    #     if(norm_method == 'Rarefaction') {
    #         microbiome_data_list <- rarefy_normalize_and_extract(sampling_depth=sample_depth,
    #                                                              filtering_quantile=low_pass_filter_threshold,
    #                                                              kraken_MRexp=data[[1]])
    #         microbiome_data_list[[7]] <- NA
    #     }
    #     else if(norm_method == 'Cumulative Sum Scaling') {
    #         microbiome_data_list <- CSS_normalize_and_extract(filtering_quantile=low_pass_filter_threshold,
    #                                                           kraken_MRexp=data[[1]])
    #         microbiome_data_list[[7]] <- NA
    #     }
    #     
    #     # Aggregation and data set creation
    #     analytic_list <- aggregate_and_filter(microbiome_data_list, data.table(metadata))
    #     
    #     # PERMANOVA table generation
    # }
    
    return(out_dt)
}









