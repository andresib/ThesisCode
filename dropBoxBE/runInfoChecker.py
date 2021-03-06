import os
import json

data = json.loads(open('/afs/cern.ch/cms/DB/conddb/test/dropbox/replay/runInfoFromLog.json','rb').read())
runInfoData = {}
emptyHltCount  = 0
emptyPromptCount  = 0
fileNameCount = 0
exportCount = 0
duplicateCount = 0
fromDuplicate = False
fromDropBoxRun = False
dropboxRunCount = 0
for filename in data:
    hlt = set()
    prompt = set()
    fromDuplicate = False
    fromDropBoxRun = False
    for (i,t,value) in data[filename]:
        value = int(value)
        if i=='dropbox_run':
            fromDropBoxRun = True
        if i=='export':
            exportCount += 1
        if i=='duplicate':
            fromDuplicate = True
        if t=='hlt' or t=='express':
            hlt.add(value)
        if t=='prompt':
            prompt.add(value)
    
    if len(hlt)>1:
        print 'HLT',data[filename],filename,hlt
        hlt = set([max(hlt)])
        print hlt
    if len(prompt)>1:
        print 'PROMPT',data[filename],filename,prompt
        prompt = set([max(prompt)])
        print prompt
    if len(hlt)==0:
        hlt.add(None)
        emptyHltCount += 1
    if len(prompt)==0:
        prompt.add(None)
        emptyPromptCount += 1
    fileNameCount += 1
    runInfoData[filename] = (hlt.pop(),prompt.pop())
    if fromDuplicate:
        duplicateCount += 1
    if fromDropBoxRun:
        dropboxRunCount += 1

with open('runInfoFromLogForReplay.json','wb') as f:
    json.dump(runInfoData,f)
print 'Files: %d empty hlt: %d empty prompt: %d. From dropbox_run %d, export: %d, from duplicate: %d.' %(fileNameCount,emptyHltCount,emptyPromptCount, dropboxRunCount, exportCount, duplicateCount )
