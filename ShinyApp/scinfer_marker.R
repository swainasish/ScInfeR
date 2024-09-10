lung<-read.csv("marker_sets/lungs.csv")
lung$tissue<-"Lungs"
pancrease<-read.csv("marker_sets/pancrease.csv",
                    row.names = 1)
pancrease$tissue<-"Pancreas"
liver<-read.csv("marker_sets/liver.csv")
liver$tissue<- "Liver"
combined_marker<- rbind(lung,pancrease,liver)
write.csv(combined_marker,
          "scinfer_combined_hs.csv")
a=data.frame(unique(liver$celltype))
write.csv(a,
          "marker_sets/delete.csv")
