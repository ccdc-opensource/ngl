
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
        expect(structure.unitcell.a).toBe(12.1359);
        expect(structure.unitcell.b).toBe(9.8351);
        expect(structure.unitcell.c).toBe(11.185);
        expect(structure.unitcell.alpha).toBe(90);
        expect(structure.unitcell.beta).toBe(90);
        expect(structure.unitcell.gamma).toBe(90);
        expect(structure.unitcell.volume).toBe(1335.016706350495);
      })
    })
  })
})
