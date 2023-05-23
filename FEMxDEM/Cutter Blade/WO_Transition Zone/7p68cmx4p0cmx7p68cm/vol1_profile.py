#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND #NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#####################################################################################
import sys; sys.path.append("..")
from profiler import *

stats = StatsWrapper('/tmp/vol1.pro')
tot = stats.getKeyByName('solveAt','demfemcoupling.py')
getForcesFromDem = stats.getKeyByName('getForcesFromDem','demfemcoupling')
applyForcesOnFem = stats.getKeyByName('applyForcesOnFem','demfemcoupling')
getDsplFromFem = stats.getKeyByName('getDsplFromFem','demfemcoupling')
applyDsplOnDem = stats.getKeyByName('applyDsplOnDem','demfemcoupling')
yadeSolve = stats.getKeyByNameAndLine('solve',193)
oofemSolve = stats.getKeyByNameAndLine('solve',135)

printer = StatsPrinter(tot.totTime,25)
printer.add('total time for solution',tot.totTime)
printer.add('getForcesFromDem',getForcesFromDem.totTime)
printer.add('applyForcesOnFem',applyForcesOnFem.totTime)
printer.add('getDsplFromFem',getDsplFromFem.totTime)
printer.add('applyDsplOnDem',applyDsplOnDem.totTime)
printer.add('yadeSolve',yadeSolve.totTime)
printer.add('oofemSolve',oofemSolve.totTime)
printer.pprint()
