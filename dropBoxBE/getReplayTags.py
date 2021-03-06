'''dropBox backend's script to get the replay tags from the original dropBox files.
'''

__author__ = 'Miguel Ojeda'
__copyright__ = 'Copyright 2012, CERN CMS'
__credits__ = ['Giacomo Govi', 'Salvatore Di Guida', 'Miguel Ojeda', 'Andreas Pfeiffer']
__license__ = 'Unknown'
__maintainer__ = 'Miguel Ojeda'
__email__ = 'mojedasa@cern.ch'


import os
import logging
import tarfile
import json

import replay
import metadata


outputFile = 'replayTags.json'


def main():
    replayTags = {}

    files = replay.getFiles()

    i = 0
    for fileName in files:
        i += 1
        logging.info('[%s/%s] %s: Getting replay tags...', i, len(files), fileName)

        tarFile = tarfile.open(os.path.join(replay.dropBoxReplayFilesFolder, fileName))

        names = tarFile.getnames()
        if len(names) != 2:
            raise Exception('%s: Invalid number of files in tar file.', fileName)

        baseFileName = names[0].rsplit('.', 1)[0]
        dbFileName = '%s.db' % baseFileName
        txtFileName = '%s.txt' % baseFileName
        if set([dbFileName, txtFileName]) != set(names):
            raise Exception('%s: Invalid file names in tar file.', fileName)
        
        newMetadata = json.loads(metadata.port(tarFile.extractfile(txtFileName).read(), fileName))
        tarFile.close()

        for tag in newMetadata['destinationTags']:
            replayTags.setdefault(newMetadata['destinationDatabase'], set([])).add(tag)
            for dependentTag in newMetadata['destinationTags'][tag]['dependencies']:
                replayTags[newMetadata['destinationDatabase']].add(dependentTag)

    # Convert the sets to a sorted list. This allows to dump the JSON and
    # also makes future changes easier to compare (since the lists are sorted,
    # when taking a diff we will not see movement of tags)
    outputDict = {}
    for destDB in replayTags:
        outputDict[destDB] = sorted(replayTags[destDB])

    logging.info('Writing output file %s...', outputFile)
    with open(outputFile, 'w') as f:
        json.dump(outputDict, f, sort_keys = True, indent = 4)


if __name__ == '__main__':
    main()

