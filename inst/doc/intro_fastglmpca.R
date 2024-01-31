## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment = "#",collapse = TRUE,results = "hold",
                      fig.align = "center",dpi = 120)

## ----load-pkgs, message=FALSE, warning=FALSE----------------------------------
library(Matrix)
library(fastglmpca)
library(ggplot2)
library(cowplot)

## ----set-seed-----------------------------------------------------------------
set.seed(1)

## ----load-data----------------------------------------------------------------
data(pbmc_facs)
dim(pbmc_facs$counts)

## ----nonzero-rate-------------------------------------------------------------
mean(pbmc_facs$counts > 0)

## ----subset-genes-1-----------------------------------------------------------
counts <- pbmc_facs$counts
n      <- nrow(counts)
rows   <- sample(n,3000)
counts <- counts[rows,]

## ----subset-genes-2-----------------------------------------------------------
dim(counts)

## ----init-glmpca--------------------------------------------------------------
fit0 <- init_glmpca_pois(counts,K = 2)

## ----fit-glmpca-1, results="hide", eval=FALSE---------------------------------
#  fit <- fit_glmpca_pois(counts,fit0 = fit0)

## ----fit-glmpca-2-------------------------------------------------------------
fit <- pbmc_facs$fit

## ----fit-glmpca-3-------------------------------------------------------------
names(fit)

## ----fitted-counts------------------------------------------------------------
fitted_counts <- with(fit,
  exp(tcrossprod(U %*% diag(d),V) +
      tcrossprod(X,B) +
      tcrossprod(W,Z)))

## ----fitted-vs-observed, fig.height=3, fig.width=3----------------------------
i <- sample(prod(dim(counts)),2e4)
pdat <- data.frame(obs    = as.matrix(counts)[i],
                   fitted = fitted_counts[i])
ggplot(pdat,aes(x = obs,y = fitted)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  labs(x = "observed count",y = "fitted count") +
  theme_cowplot(font_size = 12)

## ----plot-v, fig.height=3, fig.width=4.5--------------------------------------
celltype_colors <- c("forestgreen","dodgerblue","darkmagenta",
                     "gray","hotpink","red")
celltype <- as.character(pbmc_facs$samples$celltype)
celltype[celltype == "CD4+/CD25 T Reg" |
         celltype == "CD4+ T Helper2" |
         celltype == "CD8+/CD45RA+ Naive Cytotoxic" |
		 celltype == "CD4+/CD45RA+/CD25- Naive T" |
		 celltype == "CD4+/CD45RO+ Memory"] <- "T cell"
celltype <- factor(celltype)
pdat <- data.frame(celltype = celltype,
                   pc1 = fit$V[,1],
				   pc2 = fit$V[,2])
ggplot(pdat,aes(x = pc1,y = pc2,color = celltype)) +
  geom_point() +
  scale_color_manual(values = celltype_colors) +
  theme_cowplot(font_size = 10)

## ----session-info-------------------------------------------------------------
sessionInfo()

