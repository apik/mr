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

diffMfin.push_back(  - protHHHHH->M(0) + Mfin1(MMH,MMH,MMH,MMH,MMH) );
diffMfin.push_back(  - protHZHZZ->M(0) + Mfin1(MMH,MMZ,MMH,MMZ,MMZ) );
diffMfin.push_back(  - protHWHWW->M(0) + Mfin1(MMH,MMW,MMH,MMW,MMW) );
diffMfin.push_back(  - protHtHtt->M(0) + Mfin1(MMH,MMt,MMH,MMt,MMt) );
diffMfin.push_back(  - protHZHZZ->M(0) + Mfin1(MMZ,MMH,MMZ,MMH,MMZ) );
diffMfin.push_back(  - protZZZZH->M(0) + Mfin1(MMZ,MMZ,MMZ,MMZ,MMH) );
diffMfin.push_back(  - protZWZWW->M(0) + Mfin1(MMZ,MMW,MMZ,MMW,MMW) );
diffMfin.push_back(  - protZtZtt->M(0) + Mfin1(MMZ,MMt,MMZ,MMt,MMt) );
diffMfin.push_back(  - protHWHWW->M(0) + Mfin1(MMW,MMH,MMW,MMH,MMW) );
diffMfin.push_back(  - protZWZWW->M(0) + Mfin1(MMW,MMZ,MMW,MMZ,MMW) );
diffMfin.push_back(  - protWWWWH->M(0) + Mfin1(MMW,MMW,MMW,MMW,MMH) );
diffMfin.push_back(  - protWWWWZ->M(0) + Mfin1(MMW,MMW,MMW,MMW,MMZ) );
diffMfin.push_back(  - protWWWW0->M(0) + Mfin1(MMW,MMW,MMW,MMW,0) );
diffMfin.push_back(  - protWtWt0->M(0) + Mfin1(MMW,MMt,MMW,MMt,0) );
diffMfin.push_back(  - protHtHtt->M(0) + Mfin1(MMt,MMH,MMt,MMH,MMt) );
diffMfin.push_back(  - protZtZtt->M(0) + Mfin1(MMt,MMZ,MMt,MMZ,MMt) );
diffMfin.push_back(  - protWtWt0->M(0) + Mfin1(MMt,MMW,MMt,MMW,0) );
diffMfin.push_back(  - protttttH->M(0) + Mfin1(MMt,MMt,MMt,MMt,MMH) );
diffMfin.push_back(  - protttttZ->M(0) + Mfin1(MMt,MMt,MMt,MMt,MMZ) );
diffMfin.push_back(  - prottttt0->M(0) + Mfin1(MMt,MMt,MMt,MMt,0) );
diffVfin.push_back(  - protWWWW0->Vxzuv(0) + Vfin1(MMW,MMW,0,MMW) );
diffVfin.push_back(  - prottttt0->Vxzuv(0) + Vfin1(MMt,MMt,0,MMt) );
diffUfin.push_back(  - protHHHHH->Uzxyv(0) + Ufin1(MMH,MMH,MMH,MMH) );
diffUfin.push_back(  - protHZHZZ->Uzxyv(0) + Ufin1(MMH,MMH,MMZ,MMZ) );
diffUfin.push_back(  - protHWHWW->Uzxyv(0) + Ufin1(MMH,MMH,MMW,MMW) );
diffUfin.push_back(  - protHtHtt->Uzxyv(0) + Ufin1(MMH,MMH,MMt,MMt) );
diffUfin.push_back(  - protHZHZZ->Uuyxv(0) + Ufin1(MMZ,MMZ,MMZ,MMH) );
diffUfin.push_back(  - protZWZWW->Uzxyv(0) + Ufin1(MMZ,MMZ,MMW,MMW) );
diffUfin.push_back(  - protZtZtt->Uzxyv(0) + Ufin1(MMZ,MMZ,MMt,MMt) );
diffUfin.push_back(  - protZZ00->Uxzuv(0) + Ufin1(MMZ,MMZ,0,0) );
diffUfin.push_back(  - protHWHWW->Uuyxv(0) + Ufin1(MMW,MMW,MMW,MMH) );
diffUfin.push_back(  - protZWZWW->Uuyxv(0) + Ufin1(MMW,MMW,MMW,MMZ) );
diffUfin.push_back(  - protWtWt0->Uzxyv(0) + Ufin1(MMW,MMW,0,MMt) );
diffUfin.push_back(  - protWW00->Uxzuv(0) + Ufin1(MMW,MMW,0,0) );
diffUfin.push_back(  - protHtHtt->Uuyxv(0) + Ufin1(MMt,MMt,MMH,MMt) );
diffUfin.push_back(  - protZtZtt->Uuyxv(0) + Ufin1(MMt,MMt,MMZ,MMt) );
diffUfin.push_back(  - protWtWt0->Uuyxv(0) + Ufin1(MMt,MMt,0,MMW) );
diffTfin.push_back(  - protHZHZZ->Tyzv(0) + Tfin1(MMZ,MMZ,MMH) );
diffTfin.push_back(  - protZWZWW->Txuv(0) + Tfin1(MMZ,MMW,MMW) );
diffTfin.push_back(  - protZtZtt->Txuv(0) + Tfin1(MMZ,MMt,MMt) );
diffTfin.push_back(  - protZZ00->Txuv(0) + Tfin1(MMZ,0,0) );
diffTfin.push_back(  - protZWZWW->Tuxv(0) + Tfin1(MMW,MMZ,MMW) );
diffTfin.push_back(  - protHWHWW->Tuxv(0) + Tfin1(MMW,MMW,MMH) );
diffTfin.push_back(  - protWtWt0->Txuv(0) + Tfin1(MMW,MMt,0) );
diffTfin.push_back(  - protWW00->Txuv(0) + Tfin1(MMW,0,0) );
diffTfin.push_back(  - protZtZtt->Tyzv(0) + Tfin1(MMt,MMZ,MMt) );
diffTfin.push_back(  - protWtWt0->Tyzv(0) + Tfin1(MMt,MMW,0) );
diffTfin.push_back(  - protHtHtt->Tyzv(0) + Tfin1(MMt,MMt,MMH) );
diffSfin.push_back(  - protHtHtt->Svyz(0) + Sfin1(MMH,MMt,MMt) );
diffSfin.push_back(  - protHZHZZ->Svyz(0) + Sfin1(MMZ,MMZ,MMH) );
diffSfin.push_back(  - protZtZtt->Svyz(0) + Sfin1(MMZ,MMt,MMt) );
diffSfin.push_back(  - protHWHWW->Suxv(0) + Sfin1(MMW,MMW,MMH) );
diffSfin.push_back(  - protZWZWW->Suxv(0) + Sfin1(MMW,MMW,MMZ) );
diffSfin.push_back(  - protWWWW0->Suxv(0) + Sfin1(0,MMW,MMW) );
diffSfin.push_back(  - protWtWt0->Suxv(0) + Sfin1(0,MMW,MMt) );
diffSfin.push_back(  - prottttt0->Suxv(0) + Sfin1(0,MMt,MMt) );
