#' this function draws the manhattan plot
#' genes:get_data()$genes, including: Gene, Coordination, chr, pos, REF, HET, VAR, lethal
#' sig_gene: num_gene by 4 columns ("additive", "recessive", "dominant", "TDT")
#' type: any one of the "additive", "recessive", "dominant", "TDT"
#' log_file: character
#' test: "wG2" or "woG2"
#' bin: logical
#' effective: logical vector
#' genotype: data.frame, as in get_data()$genotype
#'
mht <- function(genes, sig_gene, type, log_file, test, bin, effective, genotype = NULL) {
  cls <- c("aquamarine4", "coral4", "cornflowerblue", "antiquewhite4", "chartreuse4",
           "chocolate4", "black", "cyan4", "darkgoldenrod", "brown4", "darkmagenta",
           "darkred", "deeppink4", "gold4", "green4", "indianred", "grey66", "lightgreen",
           "mediumorchid4", "darkorange")
  alpha <- 0.05  # alpha level
  genes$chr <- paste("chr", genes$chr, sep = "")
  chrs <- unique(genes$chr)
  bonferroni <- alpha/sum(effective)

  if (type %in% c("additive", "recessive", "dominant", "TDT")) {
    genes <- cbind(genes, sig_gene = sig_gene[, type])

    # initialize plot
    par(mar = c(5, 5, 4, 4), xpd = T)
    yhei <- 1 + round(max(c(4, -log10(alpha/sum(effective)), -log10(sig_gene[,type]))))
    plot(0:5, 0:5, ylim = c(-2, yhei + 1.5), xlim = c(-2, dim(genes)[1] + 1),
         cex.lab = 2, type = "n", ylab = "-log10(pval)", xlab = "Genomic location",
         xaxt = "n", yaxt = "n")

    # plot each chr and label
    start <- 1
    for (i in 1:length(chrs)) {
      chr <- chrs[i]
      temp <- genes[genes$chr == chr, ]
      temp <- rbind(temp[1, ], temp)  # avoid the problem when there is only one gene
      lines(start + c(0, 0:(dim(temp)[1] - 2)), -log10(temp$sig_gene), col = cls[i], lwd = 5)
      text(x = start + dim(temp)[1]/2 - 1, y = -1, lab = chr, srt = 90, cex = 2)  # label chr
      start <- start + dim(temp)[1] - 1
    }

    # axis and title
    axis(2, pos = 0, at = 0:yhei, cex.axis = 1.5)
    if (bin == T) {
      bin_lab <- "binary"
    } else {
      bin_lab <- "continuous"
    }
    title(paste("Manhattan plot:", type, "model (", test, bin_lab, ")"), cex.main = 2)

    # pval cutoff
    segments(x0 = 0, x1 = dim(genes)[1], y0 = -log10(alpha), lwd = 2, lty = 3)
    segments(x0 = 0, x1 = dim(genes)[1], y0 = -log10(bonferroni), lwd = 2, lty = 3)

    # grey lines to ease reading
    for (tmp in 1:floor(yhei)) {
      segments(x0 = 0, x1 = dim(genes)[1], y0 = tmp, lwd = 1, lty = "dashed", col = "gray95")
    }

    occupied <- round(yhei)  # for labeling significant genes
    # occupied y positions for gene name labels can only be integer/2

    for (i in 1:dim(genes)[1]) {
      # make a mark if all G3 failed for this gene
      if (!is.null(genotype) && all(genotype[i, ] == "FAILED")) {
        rect(xleft = i - 0.2, xright = i + 0.2, ybottom = -0.05, ytop = 0.05,
             col = "red", border = NA)
      }

      if (genes$sig_gene[i] >= alpha) {
        next
      }
      tmp <- -log10(genes$sig_gene[i])  # tmp is the y position for mht plot
      trial <- 1  # if a good position cannot be found, just pick one

      repeat {
        y_attempt <- tmp + runif(1, 0, yhei - tmp + 0.5)
        y_attempt <- round(y_attempt * 2)/2  # try this one
        if (y_attempt == occupied[length(occupied)]) {
          next
        }  # definitely cannot be the one just before it
        trial <- trial + 1
        if ((!y_attempt %in% occupied[max(1, length(occupied) - 4):length(occupied)]) ||
            trial > 20) {
          break
        }  # find a good position or fail to find one too many times
      }

      occupied <- c(occupied, y_attempt)  # add this gene y pos to the occupied list
      addAlpha <- function(name, alpha = 1) {
        v <- col2rgb(name)
        v <- v / 255
        rgb(v[1], v[2], v[3], alpha)
      }
      text.col <- "brown1"
      segment.col <- addAlpha("brown1", 0.5)
      text(x = i, y = y_attempt + 0.5, genes$Gene[i], pos = 3, col = text.col, cex = 1.2)  # the actual y pos is y_attemp+0.5
      segments(x0 = i, y0 = y_attempt + 0.5, y1 = tmp, col = segment.col, lwd = 2)  # vertical line
      segments(x0 = i - 0.5, x1 = i + 0.5, y0 = y_attempt + 0.5, col = segment.col, lwd = 2)  # horizontal short line
    }
  } else {
    # blank plot
    plot(0:20, 0:20, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")

    # draw legend
    for (i in 1:length(cls)) {
      segments(x0 = 2, x1 = 3, y0 = 20 - i, y1 = 20 - i, lwd = 4, col = cls[i])
      text(x = 3, y = 20 - i, lab = chrs[i], pos = 4)
    }

    # lethal gene label
    rect(xleft = 7 - 0.2, xright = 7 + 0.2, ytop = 19 + 0.2, ybottom = 19 - 0.2,
         col = "deeppink")
    rect(xleft = 7 - 0.2, xright = 7 + 0.2, ytop = 18 + 0.2, ybottom = 18 - 0.2,
         col = "darkorange")
    text(x = 7 + 0.2, y = 19, lab = "Lethal gene (count=0,pval<0.05)", pos = 4)
    text(x = 7 + 0.2, y = 18, lab = "Lethal gene (pval<0.05)", pos = 4)
  }

  return(bonferroni)
}

#' Draw manhattan plot from @param d
plot.manhattan <- function(d, main = "") {
  snapshot("plot.manhattan", "plot.Rdata")
  stopifnot(is.data.frame(d))
  stopifnot("Chrom" %in% names(d) )
  stopifnot("Position" %in% names(d) )
  stopifnot("Pval" %in% names(d) )

  d$Chrom <- as.character(d$Chrom)
  d$Position <- as.integer(as.character(d$Position))
  d$Pval <- as.numeric(d$Pval)
  if (!all(grepl("^chr", d$Chrom))) {
    d$Chrom = paste("chr", d$Chrom, sep = "")
  }

  # config
  mm10.chroms <- data.frame(chrom = paste("chr", seq(19), sep = ""),
                            length = c(195471971,
                                182113224,
                                160039680,
                                156508116,
                                151834684,
                                149736546,
                                145441459,
                                129401213,
                                124595110,
                                130694993,
                                122082543,
                                120129022,
                                120421639,
                                124902244,
                                104043685,
                                98207768,
                                94987271,
                                90702639,
                                61431566))
  Chrom.color <- rep(c("red", "blue"), 10)[1:19]
  stopifnot(all(d$Chrom %in% mm10.chroms$chrom ))

  # start plotting
  chrom.right <- cumsum(as.numeric(mm10.chroms$length))
  chrom.left <- c(0, chrom.right[-length(chrom.right)])
  chrom.mid <- 0.5 * (chrom.right + chrom.left)
  label <- paste("chr", seq(19), sep = "")
  xlim <- c(min(chrom.left), max(chrom.right))
  ylim <- c(0, max(4, -log10(d$Pval), na.rm= TRUE) * 1.1 + 4) # 4: annotate genes
  offset <- chrom.left[match(d$Chrom, mm10.chroms$chrom)] + d$Position
  col <- Chrom.color[match(d$Chrom, mm10.chroms$chrom)]
  plot(offset, -log10(d$Pval), axes = F, xlab = "Genomic location", ylab = "-log10(P)",
       xlim = xlim, ylim = ylim, col = col, main = main)
  axis(1, at = chrom.left, label = label, las = 2, lwd = 0, lwd.ticks = 1)
  axis(2)
  abline(v = c(0, chrom.right), col = "lightgray")

  alpha <- 0.05
  num.test <- sum(!is.na(d$Pval))
  bonferroni <- alpha / num.test
  abline(h = -log10(alpha), lty = "dotted", col = "lightgray")
  abline(h = -log10(bonferroni), lty = "dotted", col = "lightgray")

  d.highlight <- subset(d, Pval < alpha)
  print(d.highlight)
  for (i in seq_len(nrow(d.highlight))) {
    ## print(i)

    x <- chrom.left[match(d.highlight$Chrom[i], mm10.chroms$chrom)] + d.highlight$Position[i]
    y <- -log10(d.highlight$Pval[i])

    print(x)
    addAlpha <- function(name, alpha = 1) {
      v <- col2rgb(name)
      v <- v / 255
      rgb(v[1], v[2], v[3], alpha)
    }

    ## jitter lables
    tmp.x <- runif(1, min = -3e7, max = 3e7)
    tmp.y <- runif(1, min = 0.1, max = 2.5)
    lines(c(x, x + tmp.x), y + c(0.1, tmp.y), col = addAlpha("brown1"), lty = "dotted")
    ## adj = (0, 0.5) = (top, middle)
    text(x + tmp.x, y + 0.1 + tmp.y, labels = d.highlight$Gene[i], srt = 90, cex= 0.75, adj = c(0, 0.5))
  }
}
