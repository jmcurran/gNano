library(gNanoPkg)

resultsRoot = "/Users/jcur002/Dropbox/Code/git/gNano/systematic/"

g0.df = loadResults("g-0", resultsRoot = resultsRoot)
g0.ll = calcLogLik(g0.df, "g")

ln0.df = loadResults("ln-0", resultsRoot = resultsRoot)
ln0.ll = calcLogLik(ln0.df, "ln")

compareLL(ln0.ll, g0.ll, "ln-0", "g-0")

g1.df = loadResults("g-1", resultsRoot = resultsRoot)
g1.ll = calcLogLik(g1.df, "g")

compareLL(g0.ll, g1.ll, "g-0", "g-1")


sn1.df = loadResults("sn-1", resultsRoot = resultsRoot)
sn1.ll = calcLogLik(sn1.df, "sn")

sn2.df = loadResults("sn-2", resultsRoot = resultsRoot)
sn2.ll = calcLogLik(sn2.df, "sn")

compareLL(sn1.ll, sn2.ll, "sn-1", "sn-2")

sn3.df = loadResults("sn-3", resultsRoot = resultsRoot)
gNanoPkg:::obsExpectPlot(sn3.df, "sn")
sn3.ll = calcLogLik(sn3.df, "sn")

compareLL(sn2.ll, sn3.ll, "sn-2", "sn-3")

sn4.df = loadResults("sn-4", resultsRoot = resultsRoot)
sn4.ll = calcLogLik(sn4.df, "sn")

compareLL(sn3.ll, sn4.ll, "sn-3", "sn-4")

ln4.df = loadResults("ln-4", resultsRoot = resultsRoot)
ln4.ll = calcLogLik(ln4.df, "ln")
gNanoPkg:::obsExpectPlot(ln4.df, "ln")

compareLL(sn3.ll, ln4.ll, "sn-3", "ln-4")
