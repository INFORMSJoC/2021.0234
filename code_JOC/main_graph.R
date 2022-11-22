#result include: FDR comparisons, precision, recall and F-1 measure comparison.
source("Keyphrase_functions.R")
set.seed(541)
k<-5
result <- main.function2(k)
result$FDR_cut_mn
save(result,file = "/users/guanshenw/scratch/keyphrase/result201_500(submission_Apr27).Rdata")
Real_FDR_TR <- vector()
Real_FDR_BL <- vector()
mn_TP_Pos <- vector()
Total_pos_mn <- rep(0,6)
for (i in 1:(length(result$split_points)-1)){
  u_0_adjust <- result$u_0_adjust[(result$split_points[i]+1):result$split_points[i+1]]  
  Base_Line_adjust <- result$Base_Line_adjust[(result$split_points[i]+1):result$split_points[i+1]]  
  poster_pi_mn <- result$poster_pi_mn[(result$split_points[i]+1):result$split_points[i+1]] 
  Y <- rep(0,result$words[i])
  truth_Y <- Y
  Y[result$obs_label[[i]]] = 1
  truth <- result$truth[(result$split_key[i]+1):result$split_key[i+1]]
  truth_Y[truth] = 1
  pre.rec.auc.mn <- precision.recall.auc(Y,truth,truth_Y,poster_pi_mn,k,c=c(0.05,0.1,0.15,0.2,0.25,0.3))
  mn_TP_Pos <- rbind(mn_TP_Pos,c(pre.rec.auc.mn$FDR_tp,pre.rec.auc.mn$FDR_pos))
  Total_pos_mn <- Total_pos_mn + pre.rec.auc.mn$FDR_pos
  }
#keypharse201_500 <- sum(result$num_of_key)
#make barplot (Fig 1).
mn_0.15_article <- as.data.frame(cbind(result$TP_0.15mn_article,result$mn_comparison$check_mn_tr[,c(9,3)],result$mn_comparison$check_mn_bl[,c(9,3)],result$num_of_key,result$words))
colnames(mn_0.15_article) <- c("mean pos","mean TP","textrank pos","textrank TP","semi pos","semi TP","num of key","num of word")
mn_0.15_article$key_num_group <- ifelse(result$num_of_key<=14, 1 ,ifelse(result$num_of_key<=19, 2 ,ifelse(result$num_of_key <= 24, 3, 4)))
pro_group_quantile <- quantile(result$num_of_key/result$words,probs=c(0.25,0.5,0.75))
mn_0.15_article$key_pro_group <- ifelse(result$num_of_key/result$words<= pro_group_quantile[1], 1 ,ifelse(result$num_of_key/result$words <= pro_group_quantile[2], 2 ,ifelse(result$num_of_key/result$words <= pro_group_quantile[3], 3, 4)))

mn_0.15_key_group <- aggregate(mn_0.15_article[,c(1:8)],list(mn_0.15_article$key_num_group),sum)
mn_0.15_key_group$pos_pro <- mn_0.15_key_group$`mean pos`/mn_0.15_key_group$`num of word`
#mn_0.15_key_group <- aggregate(mn_0.15_article[mn_0.15_article[,3]>6,c(1:6)],list(mn_0.15_article[mn_0.15_article[,3]>6,]$key_num_group),sum)
mn_0.15_key_group$total_key <- aggregate(result$num_of_key,list(mn_0.15_article$key_num_group),sum)[,2]
mn_0.15_key_group$mean_precision <- mn_0.15_key_group$`mean TP`/mn_0.15_key_group$`mean pos`
mn_0.15_key_group$textrank_precision <- mn_0.15_key_group$`textrank TP`/mn_0.15_key_group$`textrank pos`
mn_0.15_key_group$semi_precision <- mn_0.15_key_group$`semi TP`/mn_0.15_key_group$`semi pos`
mn_0.15_key_group$mean_recall <- mn_0.15_key_group$`mean TP`/mn_0.15_key_group$total_key
mn_0.15_key_group$textrank_recall <- mn_0.15_key_group$`textrank TP`/mn_0.15_key_group$total_key
mn_0.15_key_group$semi_recall <- mn_0.15_key_group$`semi TP`/mn_0.15_key_group$total_key
mn_0.15_key_group$mean_f <- 2*(mn_0.15_key_group$mean_precision*mn_0.15_key_group$mean_recall)/(mn_0.15_key_group$mean_precision+mn_0.15_key_group$mean_recall)
mn_0.15_key_group$textrank_f <- 2*(mn_0.15_key_group$textrank_precision*mn_0.15_key_group$textrank_recall)/(mn_0.15_key_group$textrank_precision+mn_0.15_key_group$textrank_recall)
mn_0.15_key_group$semi_f <- 2*(mn_0.15_key_group$semi_precision*mn_0.15_key_group$semi_recall)/(mn_0.15_key_group$semi_precision+mn_0.15_key_group$semi_recall)
#mn_0.15_key_group

barplot(t(as.matrix(cbind(mn_0.15_key_group$mean_f,mn_0.15_key_group$textrank_f,mn_0.15_key_group$semi_f))),ylim=c(0,1),names.arg=c("A","B","C","D"),xlab="number of keyphrases",ylab="overall F-measure",beside = TRUE, legend=c("BSS","TR","SS"))

#make the table of keyphrase proportions (Table 4)
mn_0.15_key_pro <- aggregate(mn_0.15_article[,c(1:6,8)],list(mn_0.15_article$key_pro_group),sum)
mn_0.15_key_pro$total_key <- aggregate(result$num_of_key,list(mn_0.15_article$key_pro_group),sum)[,2]
mn_0.15_key_pro$pos_pro <- mn_0.15_key_pro$`mean pos`/mn_0.15_key_pro$`num of word`
mn_0.15_key_pro$mean_precision <- mn_0.15_key_pro$`mean TP`/mn_0.15_key_pro$`mean pos`
mn_0.15_key_pro$textrank_precision <- mn_0.15_key_pro$`textrank TP`/mn_0.15_key_pro$`textrank pos`
mn_0.15_key_pro$semi_precision <- mn_0.15_key_pro$`semi TP`/mn_0.15_key_pro$`semi pos`
mn_0.15_key_pro$mean_recall <- mn_0.15_key_pro$`mean TP`/mn_0.15_key_pro$total_key
mn_0.15_key_pro$textrank_recall <- mn_0.15_key_pro$`textrank TP`/mn_0.15_key_pro$total_key
mn_0.15_key_pro$semi_recall <- mn_0.15_key_pro$`semi TP`/mn_0.15_key_pro$total_key
mn_0.15_key_pro$mean_f <- 2*(mn_0.15_key_pro$mean_precision*mn_0.15_key_pro$mean_recall)/(mn_0.15_key_pro$mean_precision+mn_0.15_key_pro$mean_recall)
mn_0.15_key_pro$textrank_f <- 2*(mn_0.15_key_pro$textrank_precision*mn_0.15_key_pro$textrank_recall)/(mn_0.15_key_pro$textrank_precision+mn_0.15_key_pro$textrank_recall)
mn_0.15_key_pro$semi_f <- 2*(mn_0.15_key_pro$semi_precision*mn_0.15_key_pro$semi_recall)/(mn_0.15_key_pro$semi_precision+mn_0.15_key_pro$semi_recall)
mn_0.15_key_pro
