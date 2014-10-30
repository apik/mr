//
// MR - 2-loop matching and 3-loop Running, including full 2-loop EW corrections
// Copyright (C) 2014 Andrey Pikelner <pikelner@theor.jinr.ru>
//
// This file is part of MR.
//
// MR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MR.  If not, see <http://www.gnu.org/licenses/>.
//

diffMfin.push_back(  - protZHHZZ->M(0) + Mfin1(MMZ,MMH,MMH,MMZ,MMZ) );
diffMfin.push_back(  - protZZHHH->M(0) + Mfin1(MMZ,MMZ,MMH,MMH,MMH) );
diffMfin.push_back(  - protZWHWW->M(0) + Mfin1(MMZ,MMW,MMH,MMW,MMW) );
diffMfin.push_back(  - prottZtHt->M(0) + Mfin1(MMZ,MMt,MMH,MMt,MMt) );
diffMfin.push_back(  - protZWHWW->M(0) + Mfin1(MMW,MMZ,MMW,MMH,MMW) );
diffMfin.push_back(  - protWWWWH->M(0) + Mfin1(MMW,MMW,MMW,MMW,MMH) );
diffMfin.push_back(  - protWWWWZ->M(0) + Mfin1(MMW,MMW,MMW,MMW,MMZ) );
diffMfin.push_back(  - protWWWW0->M(0) + Mfin1(MMW,MMW,MMW,MMW,0) );
diffMfin.push_back(  - protWtWt0->M(0) + Mfin1(MMW,MMt,MMW,MMt,0) );
diffMfin.push_back(  - protW0W0t->M(0) + Mfin1(MMW,0,MMW,0,MMt) );
diffMfin.push_back(  - protW0W00->M(0) + Mfin1(MMW,0,MMW,0,0) );
diffMfin.push_back(  - prottZtHt->M(0) + Mfin1(MMt,MMZ,MMt,MMH,MMt) );
diffMfin.push_back(  - protWtWt0->M(0) + Mfin1(MMt,MMW,MMt,MMW,0) );
diffMfin.push_back(  - protttttH->M(0) + Mfin1(MMt,MMt,MMt,MMt,MMH) );
diffMfin.push_back(  - protttttZ->M(0) + Mfin1(MMt,MMt,MMt,MMt,MMZ) );
diffMfin.push_back(  - prottttt0->M(0) + Mfin1(MMt,MMt,MMt,MMt,0) );
diffMfin.push_back(  - prott0t0W->M(0) + Mfin1(MMt,0,MMt,0,MMW) );
diffMfin.push_back(  - protW0W0t->M(0) + Mfin1(0,MMW,0,MMW,MMt) );
diffMfin.push_back(  - protW0W00->M(0) + Mfin1(0,MMW,0,MMW,0) );
diffMfin.push_back(  - prott0t0W->M(0) + Mfin1(0,MMt,0,MMt,MMW) );
diffMfin.push_back(  - prot0000Z->M(0) + Mfin1(0,0,0,0,MMZ) );
diffMfin.push_back(  - prot0000W->M(0) + Mfin1(0,0,0,0,MMW) );
diffMfin.push_back(  - prot00000->M(0) + Mfin1(0,0,0,0,0) );
diffVfin.push_back(  - protWWWW0->Vzxyv(0) + Vfin1(MMW,MMW,0,MMW) );
diffVfin.push_back(  - prottttt0->Vzxyv(0) + Vfin1(MMt,MMt,0,MMt) );
diffUfin.push_back(  - protZHHZZ->Uzxyv(0) + Ufin1(MMH,MMZ,MMZ,MMH) );
diffUfin.push_back(  - protZWHWW->Uzxyv(0) + Ufin1(MMH,MMZ,MMW,MMW) );
diffUfin.push_back(  - prottZtHt->Uuyxv(0) + Ufin1(MMH,MMZ,MMt,MMt) );
diffUfin.push_back(  - protHZ00->Uxzuv(0) + Ufin1(MMH,MMZ,0,0) );
diffUfin.push_back(  - protZZHHH->Uyuzv(0) + Ufin1(MMZ,MMH,MMH,MMH) );
diffUfin.push_back(  - protZHHZZ->Uuyxv(0) + Ufin1(MMZ,MMH,MMZ,MMZ) );
diffUfin.push_back(  - protZWHWW->Uxzuv(0) + Ufin1(MMZ,MMH,MMW,MMW) );
diffUfin.push_back(  - prottZtHt->Uyuzv(0) + Ufin1(MMZ,MMH,MMt,MMt) );
diffUfin.push_back(  - protZWHWW->Uyuzv(0) + Ufin1(MMW,MMW,MMW,MMH) );
diffUfin.push_back(  - protZWHWW->Uuyxv(0) + Ufin1(MMW,MMW,MMW,MMZ) );
diffUfin.push_back(  - protWtWt0->Uzxyv(0) + Ufin1(MMW,MMW,0,MMt) );
diffUfin.push_back(  - protW0W00->Uzxyv(0) + Ufin1(MMW,MMW,0,0) );
diffUfin.push_back(  - protttttH->Uzxyv(0) + Ufin1(MMt,MMt,MMH,MMt) );
diffUfin.push_back(  - protttttZ->Uzxyv(0) + Ufin1(MMt,MMt,MMZ,MMt) );
diffUfin.push_back(  - protWtWt0->Uuyxv(0) + Ufin1(MMt,MMt,0,MMW) );
diffTfin.push_back(  - protZZHHH->Tuxv(0) + Tfin1(MMH,MMZ,MMH) );
diffTfin.push_back(  - protZWHWW->Tzyv(0) + Tfin1(MMH,MMW,MMW) );
diffTfin.push_back(  - prottZtHt->Tuxv(0) + Tfin1(MMH,MMt,MMt) );
diffTfin.push_back(  - protHZ00->Txuv(0) + Tfin1(MMH,0,0) );
diffTfin.push_back(  - prot0000Z->Tvxu(0) + Tfin1(MMZ,0,0) );
diffTfin.push_back(  - protWZWHW->Tzyv(0) + Tfin1(MMW,MMZ,MMW) );
diffTfin.push_back(  - protZWHWW->Tyzv(0) + Tfin1(MMW,MMW,MMH) );
diffTfin.push_back(  - protWtWt0->Tzyv(0) + Tfin1(MMW,MMt,0) );
diffTfin.push_back(  - protW0W00->Txuv(0) + Tfin1(MMW,0,0) );
diffTfin.push_back(  - prottZtHt->Tzyv(0) + Tfin1(MMt,MMZ,MMt) );
diffTfin.push_back(  - protWtWt0->Tyzv(0) + Tfin1(MMt,MMW,0) );
diffTfin.push_back(  - prottZtHt->Txuv(0) + Tfin1(MMt,MMt,MMH) );
diffSfin.push_back(  - prottZtHt->Suxv(0) + Sfin1(MMH,MMt,MMt) );
diffSfin.push_back(  - protZZHHH->Suxv(0) + Sfin1(MMZ,MMH,MMH) );
diffSfin.push_back(  - prottZtHt->Svyz(0) + Sfin1(MMZ,MMt,MMt) );
diffSfin.push_back(  - protZWHWW->Svyz(0) + Sfin1(MMW,MMW,MMH) );
diffSfin.push_back(  - protWZWHW->Svyz(0) + Sfin1(MMW,MMW,MMZ) );
diffSfin.push_back(  - protWWWW0->Suxv(0) + Sfin1(0,MMW,MMW) );
diffSfin.push_back(  - protWtWt0->Svyz(0) + Sfin1(0,MMW,MMt) );
diffSfin.push_back(  - prottttt0->Suxv(0) + Sfin1(0,MMt,MMt) );
