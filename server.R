library(shiny)
library(cummeRbund)

shinyServer(function(input, output, session) {

  output$fileinput <- renderText({
    input$action
    isolate({
      if (file.exists(paste(getwd(),"data",input$file,sep="/"))) inputfile<-paste(getwd(),"data",input$file,sep="/") 
      else inputfile<-as.character(input$file)
      inputfile
      })
    })
  
  cuff <- reactive({
    if (file.exists(paste(getwd(),"data",input$file,sep="/")))  readCufflinks(dir=paste(getwd(),"data",input$file,sep="/")  , runInfoFile="run.info",repTableFile="read_groups.info", geneFPKM = "genes.fpkm_tracking", geneDiff = "gene_exp.diff", geneCount="genes.count_tracking", geneRep="genes.read_group_tracking", isoformFPKM = "isoforms.fpkm_tracking", isoformDiff = "isoform_exp.diff", isoformCount="isoforms.count_tracking", isoformRep="isoforms.read_group_tracking", TSSFPKM = "tss_groups.fpkm_tracking", TSSDiff = "tss_group_exp.diff", TSSCount="tss_groups.count_tracking", TSSRep="tss_groups.read_group_tracking", CDSFPKM = "cds.fpkm_tracking", CDSExpDiff = "cds_exp.diff", CDSCount="cds.count_tracking", CDSRep="cds.read_group_tracking", CDSDiff = "cds.diff", promoterFile = "promoters.diff", splicingFile = "splicing.diff")
    else readCufflinks(dir= as.character(input$file) , runInfoFile="run.info",repTableFile="read_groups.info", geneFPKM = "genes.fpkm_tracking", geneDiff = "gene_exp.diff", geneCount="genes.count_tracking", geneRep="genes.read_group_tracking", isoformFPKM = "isoforms.fpkm_tracking", isoformDiff = "isoform_exp.diff", isoformCount="isoforms.count_tracking", isoformRep="isoforms.read_group_tracking", TSSFPKM = "tss_groups.fpkm_tracking", TSSDiff = "tss_group_exp.diff", TSSCount="tss_groups.count_tracking", TSSRep="tss_groups.read_group_tracking", CDSFPKM = "cds.fpkm_tracking", CDSExpDiff = "cds_exp.diff", CDSCount="cds.count_tracking", CDSRep="cds.read_group_tracking", CDSDiff = "cds.diff", promoterFile = "promoters.diff", splicingFile = "splicing.diff")
  }) 
  
  output$sampleSummary <- renderTable({
    input$action
    isolate({
      replicates(cuff())[,-1]
      })
    })
 
  output$dataSummary <- renderTable({
    input$action
    isolate({ 
      resTable <- matrix(NA, nrow=1, ncol=5)
      colnames(resTable) <- c("samples", "Gene", "isoforms", "Tss", "splicing")
      rownames(resTable) <- "No."
      resTable[1,1]  <- as.integer(dim(genes(cuff()))[2])
      resTable[1,2]  <- as.integer(dim(genes(cuff()))[1])
      resTable[1,3]  <- as.integer(dim(isoforms(cuff()))[1])
      resTable[1,4]  <- as.integer(dim(TSS(cuff()))[1])
      resTable[1,5]  <- as.integer(dim(splicing(cuff())))
      resTable
      })
  })
  
  output$mdsPlot <- renderPlot({
    genes.MDS.rep<-MDSplot(genes(cuff()),replicates=T)
    genes.MDS.rep<- genes.MDS.rep + geom_point(aes(x=M1, y=M2, fill=names,colour=names,shape=names, size=names) )
    genes.MDS.rep<- genes.MDS.rep + theme(axis.text.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.text.y=element_text(size=24, colour=1))
    genes.MDS.rep<- genes.MDS.rep + theme(axis.title.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.title.y=element_text(size=24, colour=1))
    genes.MDS.rep + theme(legend.key.width=unit(5, "cm"), legend.key.height=unit(1.5, "cm"), legend.text = element_text(colour=1, size=18, hjust=6), legend.title = element_text(colour=1, size=24), legend.title.align=1 )
    })  
  
  output$sammpleHeatmap <- renderPlot({
    if (input$var==1) { myDistHeat<-csDistHeat(genes(cuff()), replicates=T)}
    else {{ myDistHeat<-csDistHeat(genes(cuff()), replicates=F)}}
    myDistHeat <- myDistHeat + theme(axis.text.x=element_text(size=20, angle=0 ,colour=1,hjust=0.5),axis.text.y=element_text(size=20, colour=1))
    myDistHeat + theme(legend.key.width=unit(3, "cm"), legend.key.height=unit(1, "cm"), legend.text = element_text(colour=1, size=18, hjust=6), legend.title = element_text(colour=1, size=22), legend.title.align=0 )    
  })  

  output$densPlot <- renderPlot({
    if (input$var==1) { densPlot<-csDensity(genes(cuff()), replicates=T)}
    else {{ densPlot<-csDensity(genes(cuff()), replicates=F)}}
    densPlot <- densPlot + theme(axis.text.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.text.y=element_text(size=24, colour=1))
    densPlot<- densPlot + theme(axis.title.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.title.y=element_text(size=24, colour=1))
    densPlot<- densPlot + theme(legend.key.width=unit(5, "cm"), legend.key.height=unit(1.5, "cm"), legend.text = element_text(colour=1, size=18, hjust=6), legend.title = element_text(colour=1, size=24), legend.title.align=0 )
    densPlot + theme(plot.title=element_blank()) 
    
  }) 
  
  output$dispersionPlot <- renderPlot({
    dispersionPlot <- dispersionPlot(genes(cuff()))
    dispersionPlot<- dispersionPlot + theme(axis.text.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.text.y=element_text(size=24, colour=1))
    dispersionPlot<- dispersionPlot + theme(axis.title.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.title.y=element_text(size=24, colour=1))
    dispersionPlot<- dispersionPlot + theme(legend.key.width=unit(5, "cm"), legend.key.height=unit(1.5, "cm"), legend.text = element_text(colour=1, size=18, hjust=6), legend.title = element_text(colour=1, size=24), legend.title.align=0 )
    dispersionPlot + theme(strip.text=element_text(size=20)) 
  })  
   
  output$sampleboxPlot <- renderPlot({
    if (input$var==1) {boxplot <- csBoxplot(genes(cuff()), replicates=T) }
    else {boxplot <- csBoxplot(genes(cuff()), replicates=F)}
    boxplot <- boxplot + theme(axis.text.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.text.y=element_text(size=24, colour=1))
    boxplot <- boxplot + theme(axis.title.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.title.y=element_text(size=24, colour=1))
    boxplot + theme(legend.key.width=unit(5, "cm"), legend.key.height=unit(1.5, "cm"), legend.text = element_text(colour=1, size=18, hjust=6), legend.title = element_text(colour=1, size=24), legend.title.align=0 )
  }) 
  
  output$dendergramPlot <- renderPlot({
    if (input$var==0) {plot(csDendro(genes(cuff())))}
    else {plot(csDendro(genes(cuff()),replicates=T))}  
  }) 
  
  output$csScatterPlot <- renderPlot({
    if (input$var==1) { csScatterPlot<-csScatterMatrix(genes(cuff()), replicates=T)}
    else {{ csScatterPlot<-csScatterMatrix(genes(cuff()), replicates=F)}}
    csScatterPlot<- csScatterPlot + theme(axis.text.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.text.y=element_text(size=18, colour=1))
    csScatterPlot<- csScatterPlot + theme(axis.title.x=element_text(size=24, angle=0 ,colour=1,hjust=0.5),axis.title.y=element_text(size=24, colour=1))
    csScatterPlot + theme(strip.text=element_text(size=20)) 
  }) 
  
  myGenes <- reactive({ genes <- input$goiList
                        myGeneIds <- paste(unlist(strsplit(genes,split=", ")), sep=",")
                        myGenes<-getGenes(cuff(),myGeneIds)})
  myGenesNo <- reactive({ genes <- input$goiList
                        myGeneIds <- paste(unlist(strsplit(genes,split=", ")), sep=",")
                        length(myGeneIds)})
  
  output$goiHeatmap <- renderImage({
    input$goiaction
    isolate({    
      if (input$goivar==0)  hmap <- csHeatmap(myGenes(), heatscale= c(low='green',mid=1,high='red'), fullnames=T, clustering="row")
      else hmap <- csHeatmap(myGenes(), heatscale= c(low='green',mid=1,high='red'), fullnames=T, replicates=T, clustering="row")
      
             
      outfile <- tempfile(fileext='.png')
      if (myGenesNo()<5) {height=75*myGenesNo()
                          hmap <- hmap + theme(axis.ticks.x=element_blank(),axis.text.x=element_text(size=56, angle=0 ,colour=1,hjust=0.5),axis.ticks.margin=unit(0.5, "cm"),axis.text.y=element_text(size=44, colour=1), legend.text=element_text(size=36, colour=1))
                          hmap <- hmap + theme(legend.key.width=unit(5, "cm"), legend.key.height=unit(1.5, "cm"), legend.text = element_text(colour=1, size=44, hjust=3), legend.title = element_text(colour=1, size=44, hjust= 5, vjust=7) ) }
      else if (myGenesNo()>=5 & myGenesNo()<20) {height=50*myGenesNo()
                                                 hmap <- hmap + theme(axis.ticks.x=element_blank(),axis.text.x=element_text(size=56, angle=0 ,colour=1,hjust=0.5),axis.ticks.margin=unit(-1.5, "cm"),axis.text.y=element_text(size=44, colour=1), legend.text=element_text(size=36, colour=1))
                                                 hmap <- hmap + theme(legend.key.width=unit(5, "cm"), legend.key.height=unit(2, "cm"), legend.text = element_text(colour=1, size=44), legend.title = element_text(colour=1, size=44, hjust= 5, vjust=7) ) }
      else if (myGenesNo()>=20 & myGenesNo()<50) {height=25*myGenesNo()
                                                  hmap <- hmap + theme(axis.ticks.x=element_blank(),axis.text.x=element_text(size=56, angle=0 ,colour=1,hjust=0.5),axis.ticks.margin=unit(-1.5, "cm"),axis.text.y=element_text(size=44, colour=1), legend.text=element_text(size=36, colour=1))
                                                  hmap <- hmap + theme(legend.key.width=unit(7.5, "cm"), legend.key.height=unit(3, "cm"), legend.text = element_text(colour=1, size=44), legend.title = element_text(colour=1, size=44, hjust= 5, vjust=7) ) }
      else if (myGenesNo()>=50 & myGenesNo()<100) {height=20*myGenesNo()
                                                   hmap <- hmap + theme(axis.ticks.x=element_blank(),axis.text.x=element_text(size=56, angle=0 ,colour=1,hjust=0.5),axis.ticks.margin=unit(-1, "cm"),axis.text.y=element_text(size=44, colour=1), legend.text=element_text(size=36, colour=1))
                                                   hmap <- hmap + theme(legend.key.width=unit(10, "cm"), legend.key.height=unit(4, "cm"), legend.text = element_text(colour=1, size=44), legend.title = element_text(colour=1, size=44, hjust= 5, vjust=7) ) }
      else {height=20*myGenesNo()
            hmap <- hmap + theme(axis.ticks.x=element_blank(),axis.text.x=element_text(size=56, angle=0 ,colour=1,hjust=0.5),axis.ticks.margin=unit(-1.5, "cm"),axis.text.y=element_text(size=44, colour=1), legend.text=element_text(size=36, colour=1))
            hmap <- hmap + theme(legend.key.width=unit(10, "cm"), legend.key.height=unit(4, "cm"), legend.text = element_text(colour=1, size=44), legend.title = element_text(colour=1, size=44, hjust= 5, vjust=7) ) }
      
      
      if (input$goivar==0) width= 400*length(unique(replicates(cuff())$sample_name))
      else width= 200*length(replicates(cuff())$sample_name)
      
      png(outfile, width=width*3, height=height*3)
      plot(hmap+theme(strip.text.x=element_blank()) )
      dev.off()
      
      list(src = outfile,
           contentType = 'image/png',
           width=width,
           height=height,
           alt="This is alternate text")
      
    })
    }, deleteFile = TRUE) 

  
  output$goiExpression <- renderDataTable({
    input$goiaction
    isolate({
      if (input$goivar==1) { 
        myGenesTab <-repFpkm(myGenes())
        gene_short_name <- myGenes()@annotation$gene_short_name[match(myGenesTab$gene_id, myGenes()@annotation$gene_id)]
        cbind(gene_short_name, myGenesTab)
      }
      else {
        myGenesTab <-fpkm(myGenes())
        gene_short_name <- myGenes()@annotation$gene_short_name[match(myGenesTab$gene_id, myGenes()@annotation$gene_id)]
        cbind(gene_short_name, myGenesTab)
      }
    })
  })
  
  output$goiExpressionSummary <- renderTable({
    input$goiaction
    isolate({
      if (input$goivar==1) repFpkmMatrix(myGenes(), fullnames=T)
      else fpkmMatrix(myGenes(), fullnames=T)
    })
  })
 
  output$goiBoxPlot <- renderPlot({
    input$goiaction
    isolate({
      if (input$goivar==1) {
        fpkmRes <- repFpkmMatrix(myGenes(), fullnames=T)
        plotData <- as.matrix(log10(fpkmRes+1))
        par(mar=c(6,6,4,2))
        goiBoxPlot <- barplot(plotData, xlim=c(0, ((nrow(plotData)+1)*ncol(plotData)+4)), ylim=c(min(plotData), max(plotData)+2), beside=T, col=1:nrow(plotData), legend.text=T, xlab="", ylab="", axes=F, axisnames=F)
        goiBoxPlot
        mtext("samples", side=1 , line = 5, cex=2 )
        mtext("log10 FPKM", side=2 , line = 4 , cex=2 )
        axis(side=2, line=0.5, lwd=0, tck=0, cex.axis=2, col.axis=1, las=1)
        axis(side=1, line=1.5, lwd=0, tck=0, at=colMeans(goiBoxPlot),labels= colnames(plotData), cex.axis=1.5, col.axis=1, las=1)       
        axis(side=2, line=0, lwd=2, tck=-.05, col=1, labels=NA)
      }
      else {
        fpkmRes <- fpkmMatrix(myGenes(), fullnames=T)
        plotData <- as.matrix(log10(fpkmRes+1))
        par(mar=c(6,6,4,2))
        goiBoxPlot <- barplot(plotData, xlim=c(0, ((nrow(plotData)+1)*ncol(plotData)+4)), ylim=c(min(plotData), max(plotData)+2), beside=T, col=1:nrow(plotData), legend.text=T, xlab="", ylab="", axes=F, axisnames=F)
        goiBoxPlot
        mtext("samples", side=1 , line = 5, cex=2 )
        mtext("log10 FPKM", side=2 , line = 4 , cex=2 )
        axis(side=2, line=0.5, lwd=0, tck=0, cex.axis=2, col.axis=1, las=1)
        axis(side=1, line=1.5, lwd=0, tck=0, at=colMeans(goiBoxPlot),labels= colnames(plotData), cex.axis=1.5, col.axis=1, las=1)       
        axis(side=2, line=0, lwd=2, tck=-.05, col=1, labels=NA)
      }
    })
  })

  output$goiDiff <- renderDataTable({
    input$goiaction
    isolate({
      diffRes <- diffTable(myGenes())
      diffRes <- diffRes[, -length(diffRes[1,])]
      gene_short_name <- myGenes()@annotation$gene_short_name[match(diffRes$gene_id, myGenes()@annotation$gene_id)]
      cbind(gene_short_name, diffRes)
    })
  })

  output$goiFCbarplot <- renderPlot({
    input$goiaction
    isolate({
      diffDataRes <- diffData(myGenes())
      gene_short_name <- myGenes()@annotation$gene_short_name[match(diffDataRes$gene_id, myGenes()@annotation$gene_id)]
      diffDataRes <-cbind(gene_short_name, diffDataRes)
      diffDataRes <- diffDataRes[is.finite(diffDataRes$log2_fold_change),]
      par(mar=c(6,13,4,2))
      outputBarplot <- barplot(t(diffDataRes$log2_fold_change),beside=T,horiz=T, width=0.3, cex.lab=1.5, cex.axis=1.5, cex.names=1.5)
      outputBarplot
      mtext("log2 FC", side=1 , line = 4, cex=2 )
      mtext("GOI", side=2 , line = 11 , cex=2 )
      if (length(diffDataRes$gene_short_name) == length(unique(diffDataRes$gene_short_name))) {
        axis(side=2, line=-0.5, lwd=0, tck=0, at=colMeans(outputBarplot),labels=diffDataRes$gene_short_name, cex.axis=1.5, col.axis=1, las=2)
    }
      else {
        axis(side=2, line=-0.5, lwd=0, tck=0, at=colMeans(outputBarplot),labels=paste(diffDataRes$gene_short_name,diffDataRes$sample_2, diffDataRes$sample_1,sep="/"), cex.axis=1, col.axis=1, las=2)
      }
    })
  }, height = 1000, width = 800)
  
  
  output$degList <- renderDataTable({
    input$degaction
    isolate({
      degAll <- diffData(genes(cuff()))
      degAll <- degAll[degAll$status!="NOTEST",]
#      adjustp.FDR <- p.adjust(degAll$p_value, method = "fdr", n = length(degAll$p_value))
#      degAll <- as.data.frame(cbind(degAll, adjustp.FDR))      
      if (as.numeric(input$degFC)==0) {degFC <- degAll[degAll$q_value<=as.numeric(input$degFDR),]}
      else {degFC <- degAll[degAll$q_value<=as.numeric(input$degFDR) & abs(degAll$log2_fold_change) >= log2(as.numeric(input$degFC)), ]}
      gene_name <- annotation(genes(cuff()))$gene_short_name[match(degFC$gene_id, annotation(genes(cuff()))$gene_id)]
      locus <- annotation(genes(cuff()))$locus[match(degFC$gene_id, annotation(genes(cuff()))$gene_id)]
           
      degFC <- cbind(gene_name, locus, degFC)
      degFC[, -length(degFC[1,])]
    })
  })
  
  output$degSummary <- renderTable({
    input$degaction
    isolate({
      degAll <- diffData(genes(cuff()))
      if (as.numeric(input$degFC)==0) {
        degFC <- degAll[degAll$q_value<=as.numeric(input$degFDR),]
        No.DEGs<- paste(degFC$sample_2, degFC$sample_1, sep="-")   
      }
      else {
        degFC <- degAll[degAll$q_value<=as.numeric(input$degFDR) & abs(degAll$log2_fold_change) >= log2(as.numeric(input$degFC)), ]
        No.DEGs<- paste(degFC$sample_2, degFC$sample_1, sep="-")
      }
      degSummaryRes<-table(No.DEGs)
      Total<-margin.table(degSummaryRes)
      degSummaryRes<-addmargins(degSummaryRes, margin=1, FUN=sum)
      degSummaryRes
    })
  })
  
  output$cuffdiffCommand <- renderText({
    input$degaction
    isolate({
      degAll <- diffData(genes(cuff()))
      if (as.numeric(input$degFC)==0) {degFC <- degAll[degAll$q_value<=as.numeric(input$degFDR),]}
      else {degFC <- degAll[degAll$q_value<=as.numeric(input$degFDR) & abs(degAll$log2_fold_change) >= log2(as.numeric(input$degFC)), ]}
      gene_name <- annotation(genes(cuff()))$gene_short_name[match(degFC$gene_id, annotation(genes(cuff()))$gene_id)]
      if (length(gene_name)<300) {paste(gene_name, collapse=", ")}
      else {paste("more than 300 genes were identified, gene name is not presented here")}
    })    
  })
 



}
)

