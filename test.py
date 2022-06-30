## Import module
import numpy as np
import pycharmm
import pycharmm.lib as lib
import pycharmm.read as read
import pycharmm.lingo as lingo
import pycharmm.settings as settings
import pycharmm.generate as gen
import pycharmm.coor as coor
from pycharmm.cdocker import Rigid_CDOCKER

## Topology and parameter files
settings.set_bomb_level(-1)
read.rtf('"/home/trowray/acrb/Toppar/top_all36_prot.rtf"')
read.rtf('"/home/trowray/acrb/Toppar/top_all36_cgenff.rtf"', append = True)
read.prm('"/home/trowray/acrb/Toppar/par_all36m_prot.prm"', flex = True)
read.prm('"/home/trowray/acrb/Toppar/par_all36_cgenff.prm"', append = True, flex = True)
settings.set_bomb_level(0)
lingo.charmm_script('stream "/home/trowray/acrb/ligand_rtf/5qe.rtf"')

## File name and pathway
ligPDB = "/home/trowray/acrb/ligand_pdb/5qe.pdb"
ligandrtf = "/home/trowray/acrb/ligand_rtf/5qe.rtf"
confDir = "/home/trowray/acrb/ligand_conf/"
receptorPDB = "/home/trowray/acrb/pdb/1iwg.pdb"
receptorPSF = "/home/trowray/acrb/pdb/1iwg.psf"

## Build ligand and find grid box dimension
read.sequence_pdb(ligPDB)
gen.new_segment(seg_name = "LIGA")
read.pdb(ligPDB, resid = True)
coorStat = coor.stat()

xcen = coorStat.get('xave')
ycen = coorStat.get('yave')
zcen = coorStat.get('zave')

xlen = coorStat.get('xmax') - coorStat.get('xmin')
ylen = coorStat.get('ymax') - coorStat.get('ymin')
zlen = coorStat.get('zmax') - coorStat.get('zmin')

tmp = np.array([xlen, ylen, zlen])
maxlen = np.amax(tmp) + 10

## Rigid CDOCKER standard docking protocol
clusterResult, dockResult = Rigid_CDOCKER(xcen = xcen, ycen = ycen, zcen = zcen,
                                        maxlen = maxlen, ligPDB = ligPDB, receptorPDB = receptorPDB,
                                        receptorPSF = receptorPSF, confDir = confDir,
                                        flag_delete_conformer = False, flag_fast_placement = False)

print(clusterResult)
print(dockResult)
exit()