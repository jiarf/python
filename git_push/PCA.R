# === Functions Start ===
library(vegan);

PCA_analysis <- function(eset, pData, group_name, color, label=NA, filepath=NA) {
	if (is.na(filepath)) {
		filepath <- getwd();
	}
	eset <- eset[apply(abs(eset), 1, mean) > 0, ];
	pca <- prcomp(t(eset), center=TRUE, scale=FALSE);
	pca_sum <- summary(pca);
	pca_prop <- pca_sum$importance["Proportion of Variance", ];
	pca_cuml <- pca_sum$importance["Cumulative Proportion", ];
	cat("Drawing scree plot ... ");
	if (!is.na(label)) {
		filename <- paste(filepath, "/Scree_plot", label, ".pdf", sep="");
	} else {
		filename <- paste(filepath, "/Scree_plot.pdf", sep="");
	}
#	png(file=filename, height=400, width=400);
	pdf(file=filename,height=5, width=6.5);
	par(mar=c(4, 4, 1, 1), mgp=c(2.5, 0.8, 0), lwd=2);
	xlab <- "Principal component";
	ylab <- "Proportion";
	plot(0, 0, xlim=c(0.4, 10.6), ylim=c(0, 1), xlab=xlab, ylab=ylab, xaxt="n", 
		yaxt="n", xaxs="i", yaxs="i", cex.lab=1.5, col="#FFFFFF");
	axis(1, at=c(1:10), labels=c(1:10), tick=FALSE, lwd=2, cex.axis=1.5);
	axis(2, at=seq(0, 1, 0.2), labels=seq(0, 1, 0.2), lwd=2, cex.axis=1.5);
	for (i in c(1:10)) {
		rect(i - 0.4, 0, i + 0.4, pca_prop[i], col="#C5C5C5", lwd=2);
	}
	lines(c(1:10), pca_cuml[1:10], lwd=2, type="b", cex=2);
	abline(h=0.90, lty=2, lwd=2);
	rect(6.6, 0.25, 10.4, 0.15, lwd=2);
	lines(c(6.8, 7.15), c(0.22, 0.22), lwd=2);
	points(7.3, 0.22, pch=1, cex=2);
	lines(c(7.45, 7.8), c(0.22, 0.22), lwd=2);
	text(7.8, 0.22, "Cumulative", cex=1, pos=4);
	lines(c(6.8, 7.8), c(0.18, 0.18), lty=2, lwd=2);
	text(7.8, 0.18, "90% Proportion", cex=1, pos=4);
	temp <- dev.off();
	cat("Done.\n");
	
	pca_pred <- predict(pca);
	print (pca_pred)
	color_list <- c();
	for (i in c(1:ncol(eset))) {
		color_list <- c(color_list, color[[pData[pData$SampleID==colnames(eset)[i], group_name]]]);
	}
	xlab <- paste("PC 1 (",round(pca_prop[1]*100,1),'%)',sep='')
	ylab <- paste("PC 2 (",round(pca_prop[2]*100,1),'%)',sep='')
	if (!is.na(label)) {
		filename <- paste(filepath, "/PCA_plot", label, ".pdf", sep="");
	} else {
		filename <- paste(filepath, "/PCA_plot.pdf", sep="");
	}
	# prin_comp <- rda(eset, scale=TRUE);
	# ordiellipse(prin_comp, pData[, group_name], conf=0.99);
	
	cat("Drawing PCA plot ... ");
#	png(file=filename, height=400, width=480);
	pdf(file=filename, height=10, width=14);
	layout(mat=matrix(c(rep(1, 5), 2), c(1, 6)));
	par(mar=c(4, 4, 2, 2), mgp=c(2.5, 0.8, 0), lwd=2, cex=1);
	###调整ylim和xlim
	# plot(pca_pred[, 1:2], xlab=xlab, ylab=ylab, col=color_list, pch=16, 
		# cex=2, xaxt="n",  yaxt="n", cex.lab=1.5, ylim=c(-10, 16),xlim=c(-20, 20));
	#xmin=min(pca_pred[, 1])-3;xmax=max(pca_pred[, 1])+3;
	#ymin=min(pca_pred[, 2])-3;ymax=max(pca_pred[, 2])+3;
	plot(pca_pred[, 1:2], xlab=xlab, ylab=ylab, col=color_list, pch=16, 
		cex=2, xaxt="n",  yaxt="n", cex.lab=1.5,ylim=c(-10000, 20000),xlim=c(-15000, 20000));
	#print (c(ymin, ymax))
	#print (c(xmin, xmax))
	
	pData$Group <- factor(pData$Group,levels = c("case","control"))
	ordiellipse(pca_pred[, 1:2], pData$Group, kind="ehull", draw="polygon", 
		col=unique(color_list), border=unique(color_list), alpha=80, lty=2);
	axis(1, lwd=2, cex.axis=1.5);
	axis(2, lwd=2, cex.axis=1.5);
	abline(h=0, v=0, lty=2, lwd=2);
	par(mar=c(4, 0, 1, 0), cex=1);
	plot(0, 0, xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="", xaxt="n",
		yaxt="n", xaxs="i", yaxs="i", bty="n", col="#FFFFFF");
	text(0, 0.98, group_name, cex=1.5, pos=4);
	for (i in c(1:length(color))) {
		points(0.1, 0.98 - i * 0.05, col=color[[i]], pch=16, cex=2);
		text(0.15, 0.98 - i * 0.05, names(color)[i], cex=1.5, pos=4);
	}
	temp <- dev.off();
	cat("Done.\n");
}

# === Functions End ===