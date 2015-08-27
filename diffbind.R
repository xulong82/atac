library(DiffBind)

myDBA <- NULL
samples <- read.csv("samplesheet.csv",as.is=TRUE)
for(i in 1:nrow(samples)) {
  myDBA <- dba.peakset(DBA=myDBA,peaks=peaklist[[i]],
                       sampID=samples$SampleName[i], 
                       tissue=samples$Tissue[i],
                       factor=samples$Factor[i],
                       condition=samples$Condition[i],
                       treatment=samples$Treatment[i],
                       replicate=samples$Replicate[i],
                       peak.caller = "mycaller", 
                       bamReads=samples$bamReads[i],
                       control=samples$ControlID[i],
                       bamControl=samples$bamControl[i])
}
