import numpy as np
def conf2array(confs):
    import re
    pattern=re.compile(r'(\d*)([spdf])(\d*)')
    spdf={'s': 0, 'p': 1, 'd': 2, 'f': 3}
    nlocc = np.dtype([('symbol', np.unicode_, 3), ('econf', np.unicode_, 30), ('occ', np.int32, (3,))])
    data=[]
    for k,value in confs.items():
        if value is None or value.startswith('#') : continue
        vv = value.split('#')
        v=vv[0]
        m=pattern.match(v)
        n=m.group(1)
        l=spdf[m.group(2)]
        occ=m.group(3)
        if len(vv)>1 :
            econf=vv[1]
        else:
            econf=v
        if not occ : occ = 1
        data.append((k.capitalize(), econf.strip(), (n, l, occ)))
    occs = np.array(data, dtype=nlocc)
    return occs

def write_set_hubbard_n(occs, filename='set_hubbard_n.f90'):
    header="""!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
FUNCTION set_hubbard_n( psd ) RESULT( hubbard_n )
  !---------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER                      :: hubbard_n
  CHARACTER(LEN=2), INTENT(IN) :: psd
  !
  !
  SELECT CASE( TRIM(ADJUSTL(psd)) )
"""
    ender="""
     CASE DEFAULT
        !
        hubbard_n = -1
        !
        WRITE( stdout, '(/,"psd = ",A,/)' ) psd
        !
        CALL errore( 'set_hubbard_n', 'pseudopotential not yet inserted', 1 )
        !
  END SELECT
  !
  RETURN
  !
END FUNCTION set_hubbard_n
"""
    fstr='     !qepy --> auto\n'
    index=np.lexsort((occs['occ'][:,2],occs['occ'][:,1],occs['occ'][:,0]))
    occs=occs[index]
    for i in range(1,9):
        index=occs['occ'][:,0]==i
        if occs[index].size > 0 :
            fstr+='     CASE('
            for nk, s in enumerate(occs[index]['symbol']):
                if nk>0 and nk % 10 == 0: fstr += '&\n          '
                fstr += f"'{s}', "
            fstr = fstr[:-2]+')\n'
            fstr += f'        hubbard_n = {i}\n'
    fstr+='     !qepy <-- auto\n'
#             print(i, fstr)
    fstr=header+fstr+ender
    with open(filename, 'w') as fh:
        fh.write(fstr)

def write_set_hubbard_l(occs, filename='set_hubbard_l.f90'):
    header="""!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
FUNCTION set_hubbard_l( psd ) RESULT( hubbard_l )
  !---------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER                      :: hubbard_l
  CHARACTER(LEN=2), INTENT(IN) :: psd
  !
  !
  SELECT CASE( TRIM(ADJUSTL(psd)) )
"""
    ender="""
     CASE DEFAULT
        !
        hubbard_l = -1
        !
        WRITE( stdout, '(/,"psd = ",A,/)' ) psd
        !
        CALL errore( 'set_hubbard_l', 'pseudopotential not yet inserted', 1 )
        !
  END SELECT
  !
  RETURN
  !
END FUNCTION set_hubbard_l
"""
    fstr='     !qepy --> auto\n'
    index=np.lexsort((occs['occ'][:,2],occs['occ'][:,1],occs['occ'][:,0]))
    occs=occs[index]
    for i in range(1,9):
        index=occs['occ'][:,0]==i
        occs2=occs[index]
        if occs2.size > 0 :
            fstr+=f'     ! hubbard_n -> {i}\n'
            for j in range(0,4):
                index=occs2['occ'][:,1]==j
                if occs2[index].size > 0:
                    fstr+='     CASE('
                    for nk, s in enumerate(occs2[index]['symbol']):
                        fstr += f"'{s}', "
                    fstr = fstr[:-2]+')\n'
                    fstr += f'        hubbard_l = {j}\n'
    fstr+='     !qepy <-- auto\n'
    fstr=header+fstr+ender
#     print(fstr)
    with open(filename, 'w') as fh:
        fh.write(fstr)

def write_tabd(occs, filename='tabd.f90'):
    header="""!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION hubbard_occ( psd )
  !-----------------------------------------------------------------------
  !! This routine is a table (far from being complete) for the total number
  !! of localized electrons in transition metals or rare earths (PPs usually
  !! are built on non physical configurations).
  !
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=2), INTENT(IN) :: psd
  REAL(DP)                     :: hubbard_occ
  !
  SELECT CASE( TRIM(ADJUSTL(psd)) )
"""
    ender="""
     !
     !
     ! NOT INSERTED
     !
     CASE DEFAULT
        hubbard_occ = 0.d0
        call errore ('hubbard_occ', 'pseudopotential not yet inserted', 1)
     !
  END SELECT
  !
  RETURN
  !
END FUNCTION hubbard_occ
"""
    fstr='     !qepy --> auto\n'
    index=np.lexsort((occs['occ'][:,2],occs['occ'][:,1],occs['occ'][:,0]))
    occs=occs[index]
    for item in occs:
        fstr+='     CASE('
        s=item['symbol']
        fstr += f"'{s}')"
        fstr +=f'     ! {item["symbol"]}->{item["econf"]}\n'
        fstr += f'        hubbard_occ = {str(item["occ"][2])}d0\n'
    fstr+='     !qepy <-- auto\n'
    fstr=header+fstr+ender
#     print(fstr)
    with open(filename, 'w') as fh:
        fh.write(fstr)

def get_ldau_occs(filename):
    from configparser import ConfigParser
    config = ConfigParser()
    config.read(filename)
    confs=config['DEFAULT']
    occs=conf2array(confs)
    return occs

def ini2files(filename = 'qepy_econf.ini'):
    occs = get_ldau_occs(filename)
    write_set_hubbard_n(occs)
    write_set_hubbard_l(occs)
    write_tabd(occs)


if __name__ == '__main__':
    import sys
    filename = sys.argv[1] if len(sys.argv)>1 else 'qepy_econf.ini'
    ini2files(filename)
