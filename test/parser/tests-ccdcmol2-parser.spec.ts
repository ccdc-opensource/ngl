
import StringStreamer from '../../src/streamer/string-streamer'
import CCDCMol2Parser from '../../src/parser/ccdcmol2-parser'


import { join } from 'path'
import * as fs from 'fs'

describe('parser/ccdcmol2-parser', function () {
  describe('parsing', function () {
    it('basic', function () {
      var file = join(__dirname, '/../data/ABAQEB.json')
      var str = fs.readFileSync(file, 'utf-8')
      var streamer = new StringStreamer(str)
      var mol2Parser = new CCDCMol2Parser(streamer, {})
      return mol2Parser.parse().then(function (structure) {
        expect(structure.atomCount).toBe(37)
        expect(structure.bondCount).toBe(37)
      })
    })
  })
})
