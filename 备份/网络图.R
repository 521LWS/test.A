# 安装和加载RCy3包

if (!requireNamespace("RCy3", quietly = TRUE))
  BiocManager::install("RCy3",dependencies = TRUE)

library(RCy3)
cytoscapePing ()
cytoscapeVersionInfo ()

importNetworkFromFile(cell_network.sif)
importNetworkFromFile(cell_edg.txt)


help(package=RCy3)
browseVignettes("RCy3")

NET$node$id=NET$node$Key
NET$node$score=as.integer(c(1:22))
NET$network$source=NET$network$Source
NET$network$interaction=NET$network$Interaction
NET$network$target=NET$network$Target
NET$network$weight=NET$edge$Enrichment.Ratio


NET$edge$source =NET$edge$Favorable.Cell.Type
NET$edge$interaction=NET$network$Interaction
NET$edge$target =NET$edge$Unfavorable.Cell.Type
NET$edge$weight=NET$edge$Enrichment.Ratio


# From dataframes
nodes <- data.frame(id=NET$node$Key,
                    group=NET$node$Group, # categorical strings
                    score=as.integer(c(1:22)))
             
edges <- data.frame(source=NET$edge$Favorable.Cell.Type,
                    target=NET$edge$Unfavorable.Cell.Type,
                   
                    weight=as.numeric(NET$edge$Enrichment.Ratio))  # numeric
createNetworkFromDataFrames(nodes, edges, title="my first network", collection="DataFrame Example")




library(igraph)
# 创建图形对象
net=NET$network[,-2]
g <- graph_from_data_frame(net, directed = FALSE)
# 绘制网络图
plot(g, layout = layout_with_fr)
write.table(net, "_results.txt", row.names = FALSE, col.names = FALSE)






















#####转换风格可选
setVisualStyle('Marquee')
style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','group','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)
