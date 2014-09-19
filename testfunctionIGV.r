library(tractor.base)
library(XML)
SampleSheet <- cbind(
  c("Pu1","Myc","Ik_prePro","Ik_pro"),
  c("ProB","ProB","PreProB","ProB"),
  c("Pu1","Myc","Ik","Ik")
)
colnames(SampleSheet) <- c("SampleName",'Tissue',"Factor")
fileSheet <- cbind(
  c("Pu1","Myc","Ik_prePro","Ik_pro"),
  c(NA,NA,NA,NA),
  c("/Users/tcarroll/Downloads/randomTracks-2/Pu1DupMarkedNormalised.bw",
    "/Users/tcarroll/Downloads/randomTracks-2/MycDupMarkedNormalised.bw",
    "/Users/tcarroll/Downloads/randomTracks-2/Ikaros_2_preproBDupMarkedNormalised.bw",
    "/Users/tcarroll/Downloads/randomTracks-2/Ikaros_1_proBDupMarkedNormalised.bw"),
  c("/Users/tcarroll/Downloads/randomTracks/Pu1_WithInput_Input_2_proB_peaks.bed",
    "/Users/tcarroll/Downloads/randomTracks/Myc_WithInput_Input_Ch12_peaks.bed",
    "/Users/tcarroll/Downloads/randomTracks/Ikaros_2_preproB_WithInput_Input_2_proB_peaks.bed",
    "/Users/tcarroll/Downloads/randomTracks/Ikaros_1_proB_WithInput_Input_2_proB_peaks.bed")
)
colnames(fileSheet) <- c("SampleName","bam","bigwig","interval")
MakeIGVSampleMetadata(SampleSheet,fileSheet,"/Users/tcarroll/Documents")
MakeIGVSessionXML(fileSheet,"/Users/tcarroll/Documents","WeiIGVNewtest","mm9",locusName="All")

makebedtable <- function(grangesObject,name,basedirectory){

  dataTableJS <- readLines(system.file(package="tracktables","js","datatables.js"))
  jqueryJS <- readLines(system.file(package="tracktables","js","jquery.min.js"))
  dataTableCSS <- readLines(system.file(package="tracktables","js","jquery.datatables.css"))
  dataTableScroller <- readLines(system.file(package="tracktables","js","dataTables.scroller.min.js"))
  jsarray <- paste("[",paste0("[",apply(as.data.frame(grangesObject),1,function(x)paste0(shQuote(x),collapse=",")),"]",collapse=",\n"),"]")
  jsArrayForIGV <- paste0("var igvtable =",jsarray,";\n")
  jspart2 <- paste0(
    "$(document).ready(function() {
    $('#demo').html( '<table cellpadding=\"0\" cellspacing=\"0\" border=\"0\" class=\"display\" id=\"example\"></table>' );
    $('#example').dataTable( {
        deferRender:    true,
        dom:            \"frtiS\",
        scrollY:        200,
        scrollCollapse: true,

    \"data\": igvtable,\ncolumns:",
    paste0("[",paste0(
      unlist(lapply(c(colnames(as.data.frame(grangesObject))),function(x)paste0(
        c("{\"title\"",paste0(
          "\"",
          x,"\"}")
        ),collapse=":")
      )),collapse=",\n")
      ,"]")
    ,"\n","} );\n","} );\n")

  jspart1.2 <- paste0(jsArrayForIGV,jspart2)
  doc <- newXMLDoc(isHTML = T)
  html <- newXMLNode("html",parent=doc)
  head <- newXMLNode("head",parent = html)
  title <- newXMLNode("title",
                      "IGV tracktables table",
                      parent=head)
  css <- newXMLNode("style",
                    attrs=c("style type"="text/css","class"="init"),
                    paste0(dataTableCSS,collapse=""),
                    parent=head)
  jqueryjs <- newXMLNode("script",
                         attrs=c(type="text/javascript",language="javascript"),
                         paste0(jqueryJS,collapse=""),
                         parent=head)
  datatablejs <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            paste0(dataTableJS,collapse=""),
                            parent=head)
  datatableScroller <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            paste0(dataTableScroller,collapse=""),
                            parent=head)  
  jspart1.2js <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            jspart1.2,
                            parent=head)
  body <- newXMLNode("body",
                     attrs=c(class="dt-example"),
                     parent=html)
  div <- newXMLNode("div",
                    attrs=c(class="container"),
                    parent=body)
  section <- newXMLNode("section",
                        parent=div)
  h1 <- newXMLNode("h1","IGV tracktables Example",
                   parent=section)
  div2 <- newXMLNode("div",
                     attrs=c(id="demo"),
                     parent=section)
  saveXML(doc,file=file.path(basedirectory,filename),doctype="html")
  
  
}

maketracktable <- function(fileSheet,SampleSheet,filename,basedirectory,genome){

  basedirectory <- gsub("/$","",basedirectory)
  MakeIGVSampleMetadata(SampleSheet,fileSheet,basedirectory)
  
  #genome <- "mm9"
  xmlFiles <- unlist(lapply(seq(1,nrow(fileSheet)),function(x)
  MakeIGVSessionXML(fileSheet[x,,drop=F],
                    basedirectory,
                    paste0(fileSheet[x,1],"igv"),
                    genome,
                    locusName="All")
  ))

  dataTableJS <- readLines(system.file(package="tracktables","js","datatables.js"))
  jqueryJS <- readLines(system.file(package="tracktables","js","jquery.min.js"))
  dataTableCSS <- readLines(system.file(package="tracktables","js","jquery.datatables.css"))
  dataTableScroller <- readLines(system.file(package="tracktables","js","dataTables.scroller.min.js"))
  
  library(RJSONIO)
  files <- unlist(lapply(xmlFiles,function(x)relativePath(x,
                                                          gsub("//","/",file.path(basedirectory,filename))
                                                          )))
  t3mp <- "\"<a href=\\\"http://localhost:60151/load?file=\".concat(dir.concat(\"/"
  t4mp <- "\\\"\".concat(\""
  t5mp <- "</a>\")))"
  jsMat <- cbind(
    matrix(paste0("\"",as.vector(SampleSheet),"\""),ncol=ncol(SampleSheet),byrow=F),
    paste0(t3mp,files,"&merge=true",t4mp,">",SampleSheet[,1],t5mp)
  )
  setigv <- paste0("var igvtable = [",paste0(
    "[",apply(jsMat,1,function(x)paste0(
      x,collapse=","))
    ,"]\n",collapse=",")
  ,"];",sep="")

  jspart1 <- paste0("var loc = window.location.pathname;\n",
                  "var dir = loc.substring(0, loc.lastIndexOf('/'));\n",setigv,"\n")
  jspart2 <- paste0(
    "$(document).ready(function() {
    $('#demo').html( '<table cellpadding=\"0\" cellspacing=\"0\" border=\"0\" class=\"display\" id=\"example\"></table>' );
    $('#example').dataTable( {
    \"data\": igvtable,\ncolumns:",
    paste0("[",paste0(
      unlist(lapply(c(colnames(SampleSheet),"IGV"),function(x)paste0(
        c("{\"title\"",paste0(
          "\"",
          x,"\"}")
        ),collapse=":")
      )),collapse=",\n")
      ,"]")
    ,"\n","} );\n","} );\n")
  jspart1.2 <- paste0(jspart1,jspart2)
  doc <- newXMLDoc(isHTML = T)
  html <- newXMLNode("html",parent=doc)
  head <- newXMLNode("head",parent = html)
  title <- newXMLNode("title",
                      "IGV tracktables table",
                      parent=head)
  css <- newXMLNode("style",
                    attrs=c("style type"="text/css","class"="init"),
                    paste0(dataTableCSS,collapse=""),
                    parent=head)
  jqueryjs <- newXMLNode("script",
                         attrs=c(type="text/javascript",language="javascript"),
                         paste0(jqueryJS,collapse=""),
                         parent=head)
  datatablejs <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            paste0(dataTableJS,collapse=""),
                            parent=head)
  jspart1.2js <- newXMLNode("script",
                            attrs=c(type="text/javascript",language="javascript"),
                            jspart1.2,
                            parent=head)
  body <- newXMLNode("body",
                     attrs=c(class="dt-example"),
                     parent=html)
  div <- newXMLNode("div",
                    attrs=c(class="container"),
                    parent=body)
  section <- newXMLNode("section",
                        parent=div)
  h1 <- newXMLNode("h1","IGV tracktables Example",
                   parent=section)
  div2 <- newXMLNode("div",
                     attrs=c(id="demo"),
                     parent=section)
  saveXML(doc,file=file.path(basedirectory,filename),doctype="html")
  return(doc)
}
