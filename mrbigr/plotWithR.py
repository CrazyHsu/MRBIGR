drawPhyloPlotR = '''
    library(ggplot2)
    suppressMessages(library(ggtree))
    suppressMessages(library(treeio))
    options(warn=-1)
    drawPhyloPlot <- function(nwk_file, out_prefix, drop_lines, group_file, group_sep, pfmt){
        tree <- read.tree(nwk_file)
        if(drop_lines!=""){
            to_drop <- read.csv(drop_lines, sep=",", header=F)
            tree <- drop.tip(tree, as.vector(to_drop$V1))
        }

        circle_tree <- ggtree(tree, layout='fan', branch.length='none', ladderize=T, size=0.5)
        if(group_file!="") {
            group <- read.csv(group_file, header=T, sep=group_sep, row.names=1)
            names(group) <- c("subpop")
            groupInfo <- split(row.names(group), group$subpop)
            p <- groupOTU(circle_tree, groupInfo, 'subpop') + aes(color=subpop) + 
                scale_fill_brewer(palette="Set2") + 
                theme(legend.position = "right")
            ggsave(plot=p, file=paste(out_prefix, 'phylo.tree', pfmt, sep='.'), device=pfmt, height=5, width=5)
        }else{
            if(pfmt=="pdf"){
                ggsave(plot=circle_tree, file=paste(out_prefix, 'phylo.tree', pfmt, sep='.'), device=pfmt, height=5, width=5, useDingbats=FALSE)
            }else{
                ggsave(plot=circle_tree, file=paste(out_prefix, 'phylo.tree', pfmt, sep='.'), device=pfmt, height=5, width=5)
            }
        }
    }
'''

drawManhattanPlotR_2021 = '''
    suppressMessages(library(data.table))
    suppressMessages(library(dplyr))
    # library(CMplot)
    drawManhattanPlot <- function(gwas_res, thresholdi, out_prefix, sep, title, dpi, chrom, start, end, highlight_pos, highlight_text, pfmt) {
        sink('/dev/null')
        keepCol <- c('rs', 'chr', 'ps', 'p_wald')

        data <- fread(gwas_res, sep=sep, data.table=getOption("datatable.fread.datatable", FALSE), select=keepCol)
        thresholdi <- unlist(strsplit(thresholdi,","))
        if(thresholdi[1] == "None"){
            thresholdi[1] <- 1.0/dim(data)[1]
        }
        thresholdi <- as.numeric(thresholdi)
        lim <- -log10(min(data[,4])) + 2
        if(chrom!=""){
            data <- data %>% filter(chr==chrom)
            if(end > start){
                data <- data %>% filter(ps>=start, ps<=end)
            }
        }
        if(dim(data)[1]>0){
            highlight_pos <- unlist(strsplit(highlight_pos, ","))
            highlight_text <- unlist(strsplit(highlight_text, ","))
            new_highlight_pos <- NULL
            new_highlight_text <- NULL
            if(length(highlight_pos)==0){
                new_highlight_pos <- NULL
                new_highlight_text <- NULL
            }else if(length(highlight_pos)>0 && length(highlight_text)==0){
                for (i in seq(length(highlight_pos))){
                    if(highlight_pos[i] %in% data$rs){
                        new_highlight_pos <- c(new_highlight_pos, highlight_pos[i])
                        new_highlight_text <- c(new_highlight_text, highlight_text[i])
                    }
                }
                new_highlight_text <- new_highlight_pos
            }else{
                for (i in seq(length(highlight_pos))){
                    if(highlight_pos[i] %in% data$rs){
                        new_highlight_pos <- c(new_highlight_pos, highlight_pos[i])
                        new_highlight_text <- c(new_highlight_text, highlight_text[i])
                    }
                }
            }
            
            CMplot(data, plot.type='m', col=c("grey30", "grey60"), ylim=c(2, lim), threshold=thresholdi, 
                cex=c(0.5, 0.5, 0.5), signal.cex=c(0.5, 0.5, 0.5), threshold.col=c('red', 'green', 'blue'), 
                chr.den.col=NULL, amplify=TRUE, signal.pch=c(19, 19, 19), dpi=dpi, signal.col=c('red', 'green', 'blue'), 
                multracks=FALSE, LOG10=TRUE, file=pfmt, file.prefix=out_prefix, main=title, xlab="Chromosome",
                highlight=new_highlight_pos, highlight.text=new_highlight_text, highlight.text.col="black", highlight.text.cex=2)
        }
        sink()
    }
'''

drawQQplotR = '''
    suppressMessages(library(data.table))
    # library(CMplot)
    drawQQplot <- function(gwas_res, out_prefix, sep, title, dpi, threshold, pfmt) {
        sink('/dev/null')
        keepCol <- c('rs', 'chr', 'ps', 'p_wald')
        data <- fread(gwas_res, sep=sep, data.table=getOption("datatable.fread.datatable", FALSE), select=keepCol)
        CMplot(data, plot.type='q', col='grey30', conf.int.col='gray', signal.col='red', multracks=FALSE, LOG10=TRUE,
               threshold=threshold, signal.cex=0.5, file=pfmt, dpi=dpi, file.prefix=out_prefix, main=title)
        sink()
    }
'''

drawForestPlotR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(scales))
    suppressMessages(library(data.table))
    options(warn=-1)
    drawForestPlot <- function(mr_res, out_prefix, order_by_name, order_by_effect, order_by_group, group, 
        sep, height, width, pfmt){
        data <- fread(mr_res, sep=sep, data.table=getOption("datatable.fread.datatable", FALSE))
        for(mTrait in unique(data[,"mTrait"])){
            sub_data <- data[data$mTrait==mTrait,]
            for(snp in unique(sub_data[,"snp"])){
                d <- sub_data[sub_data$snp==snp,]
                if(order_by_name=="True"){
                    d <- d[order(d$pTrait),]
                }
                if(order_by_effect=="True"){
                    d <- d[order(d$effect),]
                }
                if(order_by_group=="True"){
                    groupD <- read.csv(group)
                    if(nrow(groupD)==0){
                        d$group <- 'group'
                    }else{
                        d <- merge(groupD, d, by.x ='id', by.y='pTrait')
                        d$group <- as.factor(d$group)
                        d <- d[order(d$group),]
                        names(d)[1] <- "pTrait"
                    }
                }

                d$pTrait <- factor(d$pTrait, d$pTrait)
                d$std <- sqrt((d$effect) ^ 2 / d$TMR)
                d$min <- d$effect - d$std
                d$max <- d$effect + d$std
                # effect
                if(order_by_group=="True"){
                    numColors <- length(levels(d$group))
                    myPalette <- hue_pal()(numColors)
                    p <- ggplot(data=d, aes(x=effect, y=pTrait, color=group)) +
                        geom_vline(xintercept=0, linetype="dashed", color="gray", size=1) +
                        geom_point(size=1.5, shape=15) +
                        geom_errorbarh(aes(xmin=min, xmax=max), height=.1) +
                        theme_bw() +
                        theme(legend.title = element_blank(),
                              legend.position = 'none',
                              plot.title = element_text(hjust=0.5),
                              panel.grid = element_blank(),
                              axis.line = element_line(colour='black', size=0.1),
                              axis.text.x = element_text(color='black'),
                              axis.text.y = element_text(color=myPalette[d$group]))+
                              xlab('') + ylab('') +
                              ggtitle(paste(mTrait, snp, sep='_'))
                }else{
                    p <- ggplot(data=d, aes(x=effect, y=pTrait, color=pTrait)) + 
                        geom_vline(xintercept=0, linetype="dashed", color="gray", size=1) +
                        geom_point(size=1.5, shape=15) +
                        geom_errorbarh(aes(xmin=min, xmax=max), height=.1) +
                        theme_bw() +
                        theme(legend.title = element_blank(),
                              legend.position = 'none',
                              plot.title = element_text(hjust=0.5),
                              panel.grid = element_blank(),
                              axis.line = element_line(colour='black', size=0.1),
                              axis.text.x = element_text(color='black'),
                              axis.text.y = element_text(color=hue_pal()(nrow(d))))+
                              xlab('') + ylab('') +
                              ggtitle(paste(mTrait, snp, sep='_'))
                }
                # numColors <- length(levels(d$group))
                # myPalette <- hue_pal()(numColors)
                # ggplot(data=d, aes(x=effect, y=pTrait, color=group)) +
                # p <- p + geom_vline(xintercept=0, linetype="dashed", color="gray", size=1) +
                #     geom_point(size=1.5, shape=15) +
                #     geom_errorbarh(aes(xmin=min, xmax=max), height=.1) +
                #     theme_bw() +
                #     theme(legend.title = element_blank(),
                #           legend.position = 'none',
                #           plot.title = element_text(hjust=0.5),
                #           panel.grid = element_blank(),
                #           axis.line = element_line(colour='black', size=0.1),
                #           axis.text.x = element_text(color='black'),
                #           axis.text.y = element_text(color=myPalette[d$group]))+
                #           xlab('') + ylab('') +
                #           ggtitle(paste(mTrait, snp, sep='_'))
                if(pfmt=="pdf"){
                    ggsave(file=paste(out_prefix, mTrait, snp, "MR", paste('forestplot', pfmt, sep='.'), sep='_'), 
                        plot=p, device=pfmt, height=height, width=width, useDingbats=FALSE)
                }else{
                    ggsave(file=paste(out_prefix, mTrait, snp, "MR", paste('forestplot', pfmt, sep='.'), sep='_'), 
                        plot=p, device=pfmt, height=height, width=width)
                }
            }
        }
    }
'''

drawScatterMrPlotR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(ggrepel))
    suppressMessages(library(data.table))
    options(warn=-1)
    drawScatterMrPlot <- function(mr_res, group, sig_p, order_by_pvalue, order_by_group, out_prefix, sep, height, width, pfmt){
        data <- fread(mr_res, sep=sep, data.table=getOption("datatable.fread.datatable", FALSE))
        for(mTrait in unique(data[,"mTrait"])){
            sub_data <- data[data$mTrait==mTrait,]
            for(snp in unique(sub_data[,"snp"])){
                d <- sub_data[sub_data$snp==snp,]
                if(group=="None"){
                    res <- d
                    res$group <- 'group'
                }else{
                    groupD <- read.csv(group)
                    if(nrow(groupD)==0){
                        res <- d
                        res$group <- 'group'
                    }else{
                        res <- merge(groupD, d, by.x ='id', by.y='pTrait')
                        res$group <- as.factor(res$group)
                    }
                }
                
                if(sig_p=="None"){
                    sig_p <- 0.05/nrow(res)
                }
                res$log10p <- -log10(res$pvalue)
                res[res$effect<0,'log10p'] <- -res[res$effect<0,'log10p']
                ylimMax <- max(res$log10p)+1
                ylimMin <- log10(sig_p) - 0.5
                
                if(order_by_group=="True"){
                    res <- res[order(res$group),]
                }
                if(order_by_pvalue=="True"){
                    res <- res[order(-res$log10p),]
                } 
                res$pos <- seq(1,nrow(res))
                scatter_plot <- ggplot(data=res) + aes(x=pos,y=log10p,color=group,fill=group)+
                      geom_hline(yintercept=log10(sig_p), linetype="dashed", color = "red", size=1)+
                      geom_hline(yintercept=-log10(sig_p), linetype="dashed", color = "red", size=1)+
                      geom_point(size=1.5,shape=19) + 
                      theme_bw()+
                      theme(legend.title = element_blank(),
                            plot.title = element_text(hjust = 0.5),
                            panel.grid = element_blank(),
                            axis.line = element_line(colour = 'black',size=0.1),
                            axis.text = element_text(color = 'black'),
                            axis.text.x = element_blank(),
                            axis.ticks.length.x = unit(0,'cm'))+
                      xlab('Relative position of exposures')+ylab('-log10(P-value)') 
                if(group!="None"){
                    scatter_plot <- scatter_plot + geom_label_repel(data=res[which(res$pvalue<=sig_p),], 
                        aes(label=id), color = "white", 
                        point.padding = NA, box.padding = 0.2) + 
                        guides(fill = guide_legend(override.aes = aes(color = NA)))
                }
                if(pfmt=="pdf"){
                    ggsave(paste(out_prefix, mTrait, snp, "MR", paste('scatterplot', pfmt, sep='.'), sep='_'), 
                        plot=scatter_plot, device = pfmt, height = height,width = width, useDingbats=FALSE)
                }else{
                    ggsave(paste(out_prefix, mTrait, snp, "MR", paste('scatterplot', pfmt, sep='.'), sep='_'), 
                        plot=scatter_plot, device = pfmt, height = height,width = width)
                }
                
            }
        }
    }
'''

drawScatterPsR = '''
    suppressMessages(library(ggplot2))
    options(warn=-1)
    drawScatterPs <- function(ps_file, group_file, out_prefix, ps_type, input_sep, group_sep, pfmt) {
        sink('/dev/null')
        data <- read.csv(ps_file, sep=input_sep, header=T)
        data <- data[c(1,2,3)]
        group <- read.csv(group_file, sep=group_sep)
        if(nrow(group)!=0){
                data$group <- group[match(data[, 1], group[, 1]), 2]
        }else{
            data$group <- 'one'
        }
        if(ps_type=="pca"){
            x_label <- "PC1" 
            y_label <- "PC2"
            out_file <- paste(out_prefix, 'scatter_PCA', pfmt, sep=".")
        }else{
            x_label <- "t-SNE1" 
            y_label <- "t-SNE2"
            out_file <- paste(out_prefix, 'scatter_tSNE', pfmt, sep=".")
        }
        names(data) <- c('sample','d1','d2','group')
        data <- na.omit(data)
        scatter_plot <- ggplot(data = data,aes(x=d1,y=d2,color=group))+
            geom_point(shape=19, size=1.5)+
            theme_bw() +
            theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1),
                  axis.text.y = element_text(colour = "black", size = 12, hjust =1),
                  axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                  axis.title.y = element_text(angle=90),
                  panel.border = element_rect(color = 'black', fill=NA, size = 1),
                  panel.grid = element_blank()
                 )+
            xlab(x_label)+ylab(y_label)
        if(nrow(group)==0){
            scatter_plot <- scatter_plot + theme(legend.position='none')
        }
        if(pfmt=="pdf"){
            ggsave(out_file, plot=scatter_plot, width = 4.2, height = 3, useDingbats=FALSE)
        }else{
            ggsave(out_file, plot=scatter_plot, width = 4.2, height = 3)
        }
        sink()
    }
'''

drawNetworkPlotR = '''
    options(warn = - 1)
    suppressMessages(library(network))
    suppressMessages(library(ggplot2))
    drawNetworkPlot <- function(m_ew, m_hub, fn, height, width, pfmt){
            net <- network(m_ew, matrix.type = 'edgelist', directed = F)
            m_hub <- m_hub[match(net%v%'vertex.names', m_hub$gene_id),]
            m_hub$hub_score_range <- '0-0.3'
            m_hub[m_hub$hub_score>=0.3 & m_hub$hub_score<0.8, 'hub_score_range'] <- '0.3-0.8'
            m_hub[m_hub$hub_score>=0.8, 'hub_score_range'] <- '0.8-1'
            m_hub$hub_score_range <- factor(m_hub$hub_score_range)
            m_hub$size <- 0.6
            m_hub[m_hub$hub_score>=0.8,'size'] <- 3
            col <- rev(scales::hue_pal()(3))[3-length(unique(m_hub$hub_score_range))+1:3]
            names(col) <- levels(m_hub$hub_score_range)
            net %e% "weight" <- as.vector(m_ew$weight)
            suppressMessages(GGally::ggnet2(net, edge.size='weight',  mode = 'kamadakawai', edge.color = 'gray80',
                size = m_hub$size, label = m_hub[m_hub$hub_score>=0.8, 'gene_id'],
                color.legend = 'hub_score', color = m_hub$hub_score_range, palette=col) +
                theme(legend.position = 'bottom') +
                scale_size_discrete(guide = 'none'))
            if(pfmt=="pdf"){
                suppressMessages(ggsave(fn, device=pfmt, height=height, width=width, useDingbats=FALSE))
            }else{
                suppressMessages(ggsave(fn, device=pfmt, height=height, width=width))
            }
    }
'''

drawEnrichPlotR = '''
    suppressMessages(library(clusterProfiler))
    suppressMessages(library(GO.db))
    suppressMessages(library(stringr))
    suppressMessages(library(plyr))
    suppressMessages(library(dplyr))
    suppressMessages(library(ggplot2))
    suppressMessages(library(AnnotationForge))
    drawEnrichPlot <- function(input_list, gene2go_file, bg_sep, plot_type, sigp, adjustp, out_prefix, height, width, pfmt){
        gene_list <- read.csv(input_list, header=F)
        gene2go <- read.csv(gene2go_file, sep=bg_sep, header=T)
        colnames(gene2go) <- c("GID", "GO")
        gene2go["EVIDENCE"] <- "IEA"
        gene_anno <- gene2go[, c('GID', 'GID', 'GID')]
        colnames(gene_anno) <- c("GID", "SYMBOL", "GENENAME")
        gene_anno <- gene_anno[!duplicated(gene_anno), ]

        # unlink("org.Custom.eg.db", recursive=T, force=T)
        if(!file.exists("org.Custom.eg.db")){
            sink('/dev/null')
            suppressMessages(makeOrgPackage(gene_info=gene_anno, go=gene2go, version='0.1', 
                maintainer='Custom <custom@someplace.org>', author='Custom <custom@someplace.org>', outputDir = '.', 
                tax_id='0000', genus='C', species='ustom', goTable='go', verbose=F))
            sink()
        }
        install.packages('./org.Custom.eg.db',repos=NULL,type="source",unlink=TRUE,quiet = T)
        library("org.Custom.eg.db", character.only = TRUE)

        goRes <- enrichGO(gene_list$V1, keyType = 'GID', OrgDb = "org.Custom.eg.db" , ont='ALL', 
            pvalueCutoff = 1, qvalueCutoff = 1)
        goResult <- goRes@result
        final_res <- data.frame()
        for(go_type in c('BP','MF','CC')){
            goResult_sub <- goResult[goResult$ONTOLOGY==go_type,]
            if((is.null(goResult_sub) || nrow(goResult_sub) == 0 || ncol(goResult_sub) == 0)){
              next
            }
            goRes@result <- goResult_sub
            goRes@ontology <- go_type
            tmp_sim <- simplify(goRes, cutoff=0.7)
            tmp_sim_df <- tmp_sim@result
            tmp_sim_df <- tmp_sim_df[tmp_sim_df$Count>=3,]
            if(nrow(tmp_sim_df)==0){
              next
            }
            final_res <- rbind(final_res, tmp_sim_df)
        }
        for(go_type in c('BP','MF','CC')){
            final_res_sub <- final_res[final_res$ONTOLOGY==go_type,]
            if((is.null(final_res_sub) || nrow(final_res_sub) == 0 || ncol(final_res_sub) == 0)){
                next
            }
            new_res <- new("enrichResult", result = final_res_sub, pvalueCutoff = sigp, pAdjustMethod = "BH",
                qvalueCutoff = sigp, gene = unique(unlist(lapply(goResult$geneID,function(x)unlist(strsplit(x,'/'))))),
                universe = 'Unknown', geneSets = list(), organism = 'Unknown', keytype = 'Unknown',
                ontology = "Unknown", readable = FALSE)
            if(plot_type=='barplot')
                p <- barplot(new_res) + scale_y_discrete(labels=function(x) str_wrap(x, width=60)) +
                    theme(panel.grid = element_blank())
            else if(plot_type=='dotplot')
                p <- dotplot(new_res) + scale_y_discrete(labels=function(x) str_wrap(x, width=60)) +
                    theme(panel.grid = element_blank())
            else if(plot_type=='cnetplot')
                p <- cnetplot(new_res) + scale_y_discrete(labels=function(x) str_wrap(x, width=60))
            else if(plot_type=='heatplot')
                p <- heatplot(new_res) + scale_y_discrete(labels=function(x) str_wrap(x, width=60))
            else if(plot_type=='emap')
                p <- emapplot(new_res, showCategory=20) + scale_y_discrete(labels=function(x) str_wrap(x, width=60))
            else if(plot_type=='upsetplot')
                p <- upsetplot(new_res) + scale_y_discrete(labels=function(x) str_wrap(x, width=60))
            if(pfmt=="pdf"){
                ggsave(paste(out_prefix, "enrich", plot_type, go_type, pfmt, sep='.'), plot=p, device=pfmt, height=height, width=width, useDingbats=FALSE)
            }else{
                ggsave(paste(out_prefix, "enrich", plot_type, go_type, pfmt, sep='.'), plot=p, device=pfmt, height=height, width=width)
            }
        }
    }

'''

drawHeatmapR = '''
    suppressMessages(library(data.table))
    suppressMessages(library(pheatmap))
    # options(warn=-1)
    drawHeatmap <- function(input, out_prefix, select_rows, select_cols, cluster_rows, cluster_cols, 
                        order_col_by, order_row_by, anno_row, anno_col, scale_by, show_rownames, show_colnames,
                        height, width, pfmt) {
        data <- fread(input, data.table=getOption("datatable.fread.datatable", FALSE))
        rownames(data) <- data[,1]
        data[,1] <- NULL

        if(select_rows!="None"){
            select_rows_vector <- as.vector(read.csv(select_rows, header=F, check.names=F)[,1])
            data <- data[which(rownames(data) %in% select_rows_vector), , drop=F]
        }
        if(select_cols!="None"){
            select_cols_vector <- as.vector(read.csv(select_cols, header=F, check.names=F)[,1])
            data <- data[, which(colnames(data) %in% select_cols_vector), drop=F]
        }

        if(order_col_by!="None"){
            order_col_by <- read.csv(order_col_by, check.names=F)
            order_col_by <- as.vector(order_col_by[,1])
            data <- data[, order_col_by, drop=F]
        }
        if(order_row_by!="None"){
            order_row_by <- read.csv(order_row_by, check.names=F)
            order_row_by <- as.vector(order_row_by[,1])
            data <- data[order_row_by, , drop=F]
        }
        if(anno_col!="None"){
            annotation_col <- read.csv(anno_col, check.names=F)
            rownames(annotation_col) <- annotation_col[,1]
            annotation_col <- annotation_col[,2:dim(annotation_col)[2], drop=F]
        }else{
            annotation_col <- NA
        }
        if(anno_row!="None"){
            annotation_row <- read.csv(anno_row, check.names=F)
            # print(annotation_row)
            rownames(annotation_row) <- annotation_row[,1]
            # print(head(annotation_row))
            annotation_row <- annotation_row[,2:dim(annotation_row)[2], drop=F]
        }else{
            annotation_row <- NA
        }
        cluster_rows <- as.logical(cluster_rows)
        cluster_cols <- as.logical(cluster_cols)
        show_rownames <- as.logical(show_rownames)
        show_colnames <- as.logical(show_colnames)
        if(show_rownames & dim(data)[1] <= 50){
            show_rownames <- T
        }else{
            show_rownames <- F
        }
        if(show_colnames & dim(data)[2] <= 50){
            show_colnames <- T
        }else{
            show_colnames <- F
        }
        out_file <- paste(out_prefix, paste("heatmap", pfmt, sep="."), sep="_")
        data_rownames <- rownames(data)
        data <- data.frame(lapply(data, function(x) {gsub("-", 0, x)}), check.names=F)
        data <- data.frame(lapply(data, as.numeric), check.names=F)
        rownames(data) <- data_rownames
        pheatmap(data, filename=out_file, scale=scale_by, cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                 show_rownames=show_rownames, show_colnames=show_colnames, annotation_row=annotation_row, 
                 annotation_col=annotation_col, angle_col=315,
                 annotation_names_row=F, annotation_names_col=F, height=height, width=width)
    }
'''

drawBoxplotR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(data.table))
    suppressMessages(library(ggpubr))
    suppressMessages(library(rstatix))
    suppressMessages(library(R.utils))
    options(warn=-1)
    removeHeterozygous <- function(df, sep="", col=1){
        g2 <- as.matrix(df[,col,drop=F])
        n <- nrow(g2)
        m <- ncol(g2)
        a1<-matrix(sapply(strsplit(g2,sep),"[",1),nrow = n,ncol=m,byrow = F)
        a2<-matrix(sapply(strsplit(g2,sep),"[",2),nrow = n,ncol=m,byrow = F)
        colnames(a1) <- colnames(g2)
        rownames(a1) <- rownames(g2)
        H<-matrix(as.numeric(!a1==a2),nrow = n,ncol = m,byrow = F)
        H <- as.data.frame(H)
        # colnames(H) <- colnames(df)
        # rownames(H) <- rownames(df)
        return(df[which(H[,1]==0), ,drop=F])
    }

    drawBoxplot <- function(input_mat_file, snp_file, group_file, selected_genes, selected_snps, selected_strains, to_scale, plot_type, out_prefix, pfmt){
        input_mat <- fread(input_mat_file, data.table=getOption("datatable.fread.datatable", FALSE))
        snp_mat <- fread(snp_file, data.table=getOption("datatable.fread.datatable", FALSE))
        group_mat <- fread(group_file, data.table=getOption("datatable.fread.datatable", FALSE))

        rownames(input_mat) <- input_mat[,1]
        input_mat[,1] <- NULL
        rownames(snp_mat) <- snp_mat[,1]
        snp_mat[,1] <- NULL
        rownames(group_mat) <- group_mat[,1]
        group_mat[,1] <- NULL

        if(selected_genes!="None"){
            selected_genes <- trimws(unlist(strsplit(selected_genes, ",")))
        }else{
            selected_genes <- c(colnames(input_mat)[1])
        }
        if(selected_snps!="None"){
            selected_snps <- trimws(unlist(strsplit(selected_snps, ",")))
        }else{
            selected_snps <- c(colnames(snp_mat)[1])
        }
        if(selected_strains!="None"){
            selected_strains <- trimws(unlist(strsplit(selected_strains, ",")))
        }else{
            selected_strains <- rownames(input_mat)
        }

        sub_input_mat <- input_mat[selected_strains, selected_genes, drop=F]
        sub_rownames <- rownames(sub_input_mat)
        sub_input_mat <- data.frame(lapply(sub_input_mat, function(x) {gsub("-", "0", x)}))
        rownames(sub_input_mat) <- sub_rownames
        sub_snp_mat <- snp_mat[selected_strains, selected_snps, drop=F]
        sub_group_mat <- group_mat[selected_strains, , drop=F]
        colnames(sub_group_mat) <- c("Subpop")

        for(snp in selected_snps){
            genotype <- sub_snp_mat[, snp, drop=F]
            # genotype <- removeHeterozygous(genotype, col=1)
            # print(head(genotype))
            colnames(genotype) <- c("Genotype")
            if(plot_type=="single_gene"){
                for(gene in selected_genes){
                    single_gene_mat <- sub_input_mat[, gene, drop=F]
                    merged_exp <- merge(x=genotype, y=single_gene_mat, by=0)
                    merged_exp <- na.omit(merged_exp)
                    colnames(merged_exp) <- c("RIL", "Genotype", "Expression")
                    merged_exp$Expression <- as.numeric(as.vector(merged_exp$Expression))
                    if(to_scale=="True"){
                        merged_exp$Expression <- (merged_exp$Expression-mean(merged_exp$Expression))/sd(merged_exp$Expression)
                    }

                    merged_exp$Genotype <- as.factor(merged_exp$Genotype)
                    p <- ggboxplot(data=merged_exp, x="Genotype", y="Expression", fill="Genotype", width=0.5, bxp.errorbar=T) +
                        theme_bw() +
                        theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                            axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                            axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                            axis.title.y = element_text(angle=90),
                            panel.border = element_rect(color = 'black', fill=NA, size = 1),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                        ) + xlab('')

                    my_comparisons <- unlist(lapply(2, combn, x = as.vector(unique(merged_exp$Genotype)), simplify = FALSE), recursive = F)
                    p <- p + stat_compare_means(method="t.test", comparisons = my_comparisons, label="p", tip.length=0.02) +
                        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
                    if(pfmt=="pdf"){
                        ggsave(paste(out_prefix, gene, paste('boxplot', pfmt, sep='.'), sep='_'), plot=p, width = 3, height = 3, useDingbats=FALSE)
                    }else{
                        ggsave(paste(out_prefix, gene, paste('boxplot', pfmt, sep='.'), sep='_'), plot=p, width = 3, height = 3)
                    }
                    
                }
            }else if(plot_type=="multi_genes") {
                gene_n <- length(selected_genes)
                if(gene_n==1){
                    width <- 3
                }else if(gene_n==2){
                    width <- 3.5
                }else{
                    width <- 3.5 + 0.5 * (gene_n-2)
                }
                if(width >= 8){
                    width <- 8
                }

                merged_exp <- merge(x=genotype, y=sub_input_mat, by=0)
                merged_exp <- na.omit(merged_exp)
                colnames(merged_exp)[1:2] <- c("RIL", "Genotype")
                if(to_scale=="True"){
                    tmp_exp <- merged_exp[,3:dim(merged_exp)[2]]
                    merged_exp[,3:dim(merged_exp)[2]] <- sapply(tmp_exp, function(tmp_exp) (tmp_exp-mean(tmp_exp))/sd(tmp_exp))
                }

                setDT(merged_exp)
                melt_merged_exp <- melt(merged_exp, id.vars = c("RIL", "Genotype"), variable.name = "Gene", value.name = "Expression")                
                melt_merged_exp$Gene <- as.factor(melt_merged_exp$Gene)
                melt_merged_exp$Genotype <- as.factor(melt_merged_exp$Genotype)
                melt_merged_exp$Expression <- as.numeric(as.vector(melt_merged_exp$Expression))
                p <- ggboxplot(data=melt_merged_exp, x="Gene", y="Expression", fill="Genotype", width=0.5, bxp.errorbar=T) + 
                    theme_bw() +
                    theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                        axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                        axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                        axis.title.y = element_text(angle=90),
                        panel.border = element_rect(color = 'black', fill=NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                    ) + xlab('')
                stat.test <- melt_merged_exp %>% group_by(Gene) %>% t_test(Expression ~ Genotype)
                stat.test <- stat.test %>% add_xy_position(x = "Gene", dodge = 0.8)
                p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02) + 
                    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
                if(pfmt=="pdf"){
                    ggsave(paste(out_prefix, 'multiGenes', snp, paste('boxplot', pfmt, sep='.'), sep='_'), plot=p, width = width, height = 3, useDingbats=FALSE)
                }else{
                    ggsave(paste(out_prefix, 'multiGenes', snp, paste('boxplot', pfmt, sep='.'), sep='_'), plot=p, width = width, height = 3)
                }
            }else if(plot_type=="single_pop") {
                for(gene in selected_genes){
                    single_gene_mat <- sub_input_mat[, gene, drop=F]
                    colnames(single_gene_mat) <- c("Expression")
                    sub_group_mat_tmp <- sub_group_mat
                    genotype_tmp <- genotype
                    single_gene_mat_tmp <- single_gene_mat

                    sub_group_mat_tmp["RIL"] <- rownames(sub_group_mat_tmp)
                    genotype_tmp["RIL"] <- rownames(genotype_tmp)
                    single_gene_mat_tmp["RIL"] <- rownames(single_gene_mat_tmp)
                    df_lst <- list(sub_group_mat_tmp, genotype_tmp, single_gene_mat_tmp)

                    merged_exp <- Reduce(function(x,y) merge(x,y,by="RIL"), df_lst)
                    merged_exp <- na.omit(merged_exp)
                    merged_exp$Expression <- as.numeric(as.vector(merged_exp$Expression))
                    if(to_scale=="True"){
                        merged_exp$Expression <- (merged_exp$Expression-mean(merged_exp$Expression))/sd(merged_exp$Expression)
                    }
                    merged_exp = removeHeterozygous(merged_exp, col="Genotype")
                    merged_exp$Subpop <- as.factor(merged_exp$Subpop)
                    merged_exp$Genotype <- as.factor(merged_exp$Genotype)
                    p <- ggboxplot(data=merged_exp, x="Subpop", y="Expression", fill="Genotype", width=0.5, 
                        bxp.errorbar=T, outlier.shape=20)+
                        theme_bw() +
                        theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                            axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                            axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                            axis.title.y = element_text(angle=90),
                            panel.border = element_rect(color = 'black', fill=NA, size = 1),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                        ) + xlab('') + coord_cartesian(ylim = c(0, NA)) 
                    stat.test <- merged_exp %>% group_by(Subpop) %>% t_test(Expression ~ Genotype)
                    stat.test <- stat.test %>% add_xy_position(x = "Subpop", dodge = 0.8)
                    p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02, size=3) + 
                        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
                    if(pfmt=="pdf"){
                        ggsave(paste(out_prefix, gene, snp, paste('pop.boxplot', pfmt, sep='.'), sep='_'), plot=p, width = 3, height = 3, useDingbats=FALSE)
                    }else{
                        ggsave(paste(out_prefix, gene, snp, paste('pop.boxplot', pfmt, sep='.'), sep='_'), plot=p, width = 3, height = 3)
                    }
                }
            }
        }
    }
'''

drawGroupPlotR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(data.table))
    suppressMessages(library(ggpubr))
    suppressMessages(library(dplyr))
    suppressMessages(library(rstatix))
    # suppressMessages(library(R.utils))
    options(warn=-1)
    drawGroupPlot <- function(input_file, sep, group_field, group_order_file, value_field, fill_field, add_pvalue, out_prefix, height, width, filter_by_count, pfmt){
        sink('/dev/null')
        input_mat <- fread(input_file, data.table=getOption("datatable.fread.datatable", FALSE))
        input_mat[, group_field] <- as.factor(input_mat[, group_field])
        if(fill_field!="None"){
            input_mat[, fill_field] <- as.factor(input_mat[, fill_field])
        }else{
            fill_field <- group_field
        }
        if(filter_by_count>0) {
            input_mat_n <- input_mat %>% group_by(!!as.name(group_field)) %>% summarise(n = n()) %>% filter(n>=filter_by_count)
            input_mat <- input_mat %>% filter(!!as.name(group_field) %in% input_mat_n[[group_field]])
        }
        if(group_order_file!="None"){
            group_order <- read.csv(group_order_file, header=FALSE)
            group_order$V1 <- as.vector(group_order$V1)
            input_mat <- input_mat %>% filter(!!as.name(group_field) %in% group_order$V1)
            input_mat[[group_field]] <- factor(input_mat[[group_field]], levels=group_order$V1)
        }
        input_mat <- input_mat %>% select(c(group_field, value_field, fill_field))
        colnames(input_mat) <- c("group", "value", "fill")
        p <- ggboxplot(data=input_mat, x="group", y="value", fill="fill", width=0.7, 
                        bxp.errorbar=T, outlier.shape=20)+
                        theme_bw() +
                        theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1, hjust = 1, angle=45),
                            axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                            axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                            axis.title.y = element_text(angle=90),
                            panel.border = element_rect(color = 'black', fill=NA, size = 1),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                        ) + 
                        xlab('')
        if(add_pvalue=="True"){
            stat.test <- input_mat %>% group_by(group) %>% t_test(value ~ fill)
            stat.test <- stat.test %>% add_xy_position(x = "group", dodge = 0.8)
            p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02)
        }
        p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
        if(pfmt=="pdf"){
            ggsave(paste(out_prefix, paste('group_box', pfmt, sep='.'), sep='_'), plot=p, width = width, height = height, useDingbats=FALSE)
        }else{
            ggsave(paste(out_prefix, paste('group_box', pfmt, sep='.'), sep='_'), plot=p, width = width, height = height)
        }
        sink()
    }
'''

drawAlleleFreqPlotR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(dplyr))
    suppressMessages(library(data.table))
    options(warn=-1)
    removeHeterozygous <- function(df, sep="", col=1){
        g2 <- as.matrix(df[,col,drop=F])
        n <- nrow(g2)
        m <- ncol(g2)
        a1<-matrix(sapply(strsplit(g2,sep),"[",1),nrow = n,ncol=m,byrow = F)
        a2<-matrix(sapply(strsplit(g2,sep),"[",2),nrow = n,ncol=m,byrow = F)
        colnames(a1) <- colnames(g2)
        rownames(a1) <- rownames(g2)
        H<-matrix(as.numeric(!a1==a2),nrow = n,ncol = m,byrow = F)
        H <- as.data.frame(H)
        # colnames(H) <- colnames(df)
        # rownames(H) <- rownames(df)
        return(df[which(H[,1]==0), ,drop=F])
    }

    drawAlleleFreqPlot <- function(input_file, snps, out_prefix, height, width, pfmt){
        sink('/dev/null')
        all_genotype <- fread(input_file, data.table=getOption("datatable.fread.datatable", FALSE), header=TRUE)
        colnames(all_genotype)[1] <- "line"
        all_genotype$group <- as.factor(all_genotype$group)
        all_snps <- unlist(strsplit(snps, ","))

        for(snp in all_snps){
            sub_genotype <- all_genotype[, c("line", "group", snp)]
            print(head(sub_genotype))
            colnames(sub_genotype) <- c("line", "group", "Genotype")
            sub_genotype$Genotype <- as.factor(sub_genotype$Genotype)
            sub_genotype <- removeHeterozygous(sub_genotype, col=3)
            sub_genotype <- sub_genotype %>% group_by(group, Genotype) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))

            p <- ggplot(sub_genotype, aes(x=group, fill = Genotype, group = Genotype)) + 
                geom_bar(aes(y=freq), stat="identity", position = "dodge") +
                theme_bw() +
                theme(axis.title.x=element_blank(),
                    axis.text.x = element_text(colour = "black", size = 12, vjust =1),
                    axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                    axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                    axis.title.y = element_text(angle=90),
                    panel.border = element_rect(color = 'black', fill=NA, size = 1),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                ) +
                scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
                ylab("Allele frequency")
            if(pfmt == "pdf"){
                ggsave(paste(out_prefix, snp,"allele_freq", pfmt, sep="."), plot=p, width = width, height = height, useDingbats=FALSE)
            }else{
                ggsave(paste(out_prefix, snp,"allele_freq", pfmt, sep="."), plot=p, width = width, height = height)
            }
        }
        sink()
    }

'''

calPcaR = '''
    suppressMessages(library(data.table))
    suppressMessages(library(dplyr))
    suppressMessages(library(FactoMineR))
    suppressMessages(library(factoextra))
    options(warn=-1)
    calPCA <- function(phe_file, group_file, selected_strains, input_sep, group_sep, group_min_n, out_prefix, plot, pfmt="pdf"){
        sink('/dev/null')

        phe_df <- fread(phe_file, data.table=getOption("datatable.fread.datatable", FALSE), check.names=F, sep=input_sep, header=TRUE)
        group_df <- fread(group_file, data.table=getOption("datatable.fread.datatable", FALSE), check.names=F, sep=group_sep, header=TRUE)
        colnames(phe_df)[1] <- "ID"
        colnames(group_df) <- c("line", "group")

        if(selected_strains!="None"){
            selected_strains <- trimws(unlist(strsplit(selected_strains, ",")))
        }else{
            selected_strains <- colnames(phe_df)
        }
        phe_df <- phe_df %>% filter(ID %in% selected_strains)

        if(group_min_n>0){
            group_df_by_n <- group_df %>% group_by(group) %>% summarise(n = n()) %>% filter(n>=group_min_n)            
            group_df <- group_df %>% filter(group %in% group_df_by_n$group)
        }
        groups <- unique(group_df$group)

        final_res <- list()
        for (var in groups){
            selected_phe <- as.vector(group_df[which(group_df$group==var),]$line)
            sub_phe_df <- phe_df[,c("ID",selected_phe), drop=F]
            rownames(sub_phe_df) <- sub_phe_df[, 1]
            sub_phe_df[, 1] <- NULL

            res.pca <- PCA(sub_phe_df, graph = FALSE)
            pc1 <- get_pca_ind(res.pca)$coord[,1,drop=F]
            colnames(pc1) <- c(var)
            final_res[[var]] <- pc1
        }
        pc_res <- Reduce(merge, lapply(final_res, function(x) data.frame(x, ID = row.names(x), check.names=FALSE)))
        write.csv(pc_res, file=paste0(out_prefix,".PCA.csv"), quote=FALSE, row.names=FALSE)
        sink()
    }

'''