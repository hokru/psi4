/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psapt.h"
#include "psi4/physconst.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libqt/qt.h"
#include <cmath>


namespace psi {
namespace sapt {

PSAPT::PSAPT(SharedWavefunction Dimer, SharedWavefunction MonomerA, SharedWavefunction MonomerB, Options& options,
                 std::shared_ptr<PSIO> psio)
    : SAPT2p3(Dimer, MonomerA, MonomerB, options, psio)
      // e_elst13_(0.0),
      // e_ind30_(0.0),
      // e_exch_ind30_(0.0),
      // e_ind30r_(0.0),
      // e_exch_ind30r_(0.0),
      // e_ind_disp30_(0.0),
      // e_exch_ind_disp30_(0.0),
      // e_disp30_(0.0),
      // e_exch_disp30_(0.0),
      // e_sapt2pp3_(0.0),
      // e_sapt2p3_(0.0),
      // e_sapt2pp3_ccd_(0.0),
//       e_sapt2p3_ccd_(0.0) {
    {third_order_ = options_.get_bool("DO_THIRD_ORDER");}
// {}

PSAPT::~PSAPT() {}

double PSAPT::compute_energy() {
    print_header();
    
    timer_on("DF Integrals       ");
    df_integrals();
    timer_off("DF Integrals       ");
    timer_on("Omega Integrals    ");
    w_integrals();
    timer_off("Omega Integrals    ");
    timer_on("Amplitudes         ");
    SAPT2p3::amplitudes();
    timer_off("Amplitudes         ");
    timer_on("Elst10             ");
    elst10();
    timer_off("Elst10             ");
    timer_on("Exch10 S^2         ");
    exch10_s2();
    timer_off("Exch10 S^2         ");
    timer_on("Exch10             ");
    exch10();
    timer_off("Exch10             ");
    timer_on("Ind20,r            ");
    ind20r();
    timer_off("Ind20,r            ");
    timer_on("Exch-Ind20,r       ");
    exch_ind20r();
    timer_off("Exch-Ind20,r       ");
    // timer_on("Disp20             ");
    // disp20();
    // timer_off("Disp20             ");
    timer_on("Exch-Disp20        ");
    exch_disp20();
    timer_off("Exch-Disp20        ");
    timer_on("Elst12             ");
    elst12();
    timer_off("Elst12             ");
    timer_on("Exch11             ");
    exch11();
    timer_off("Exch11             ");
    timer_on("Exch12             ");
    exch12();
    timer_off("Exch12             ");
    timer_on("Ind22              ");
    ind22();
    timer_off("Ind22              ");
    // timer_on("Disp21             ");
    // disp21();
    // timer_off("Disp21             ");

    // if (mbpt_disp_) {
    //     timer_on("Disp22 (SDQ)       ");
    //     disp22sdq();
    //     timer_off("Disp22 (SDQ)       ");
    //     timer_on("Disp22 (T)         ");
    //     disp22t();
    //     timer_off("Disp22 (T)         ");
    // }

    // if (ccd_disp_) {
    //     timer_on("Disp2(CCD)         ");
    //     disp2ccd();
    //     timer_off("Disp2(CCD)         ");
    //     timer_on("Disp22 (T) (CCD)   ");
    //     disp22tccd();
    //     timer_off("Disp22 (T) (CCD)   ");
    // }

    timer_on("Elst13             ");
    elst13();
    timer_off("Elst13             ");
    // timer_on("Disp30             ");
    // disp30();
    // timer_off("Disp30             ");
    if (third_order_) {
        // timer_on("ExchDisp30         ");
        // exch_disp30();
        // timer_off("ExchDisp30         ");
        timer_on("Ind30              ");
        ind30();
        timer_off("Ind30              ");
        timer_on("Ind30,r            ");
        ind30r();
        timer_off("Ind30,r            ");
        timer_on("Exch-Ind30         ");
        exch_ind30();
        timer_off("Exch-Ind30         ");
        // timer_on("IndDisp30          ");
        // ind_disp30();
        // timer_off("IndDisp30          ");
        // timer_on("ExchIndDisp30      ");
        // exch_ind_disp30();
        // timer_off("ExchIndDisp30      ");
    }

    print_results();

    return (e_saptd4_);
}

void PSAPT::print_header() {

    outfile->Printf("       SAPT-D4   \n");
    outfile->Printf("       H. Kruse, E. Hohenstein \n");
    outfile->Printf("       summer 2019\n");
    outfile->Printf("\n");
    outfile->Printf("      Orbital Information\n");
    outfile->Printf("  --------------------------\n");
    if (nsoA_ != nso_ || nsoB_ != nso_) {
        outfile->Printf("    NSO        = %9d\n", nso_);
        outfile->Printf("    NSO A      = %9zu\n", nsoA_);
        outfile->Printf("    NSO B      = %9zu\n", nsoB_);
        outfile->Printf("    NMO        = %9d\n", nmo_);
        outfile->Printf("    NMO A      = %9zu\n", nmoA_);
        outfile->Printf("    NMO B      = %9zu\n", nmoB_);
    } else {
        outfile->Printf("    NSO        = %9d\n", nso_);
        outfile->Printf("    NMO        = %9d\n", nmo_);
    }
    outfile->Printf("    NRI        = %9zu\n", ndf_);
    outfile->Printf("    NOCC A     = %9zu\n", noccA_);
    outfile->Printf("    NOCC B     = %9zu\n", noccB_);
    outfile->Printf("    FOCC A     = %9zu\n", foccA_);
    outfile->Printf("    FOCC B     = %9zu\n", foccB_);
    outfile->Printf("    NVIR A     = %9zu\n", nvirA_);
    outfile->Printf("    NVIR B     = %9zu\n", nvirB_);
    outfile->Printf("\n");

    auto mem = (long int)memory_;
    mem /= 8L;
    long int occ = noccA_;
    if (noccB_ > noccA_) occ = noccB_;
    long int vir = nvirA_;
    if (nvirB_ > nvirA_) vir = nvirB_;
    long int ovov = occ * occ * vir * vir;
    long int vvnri = vir * vir * ndf_;
    double memory = 8.0 * (vvnri + ovov * 3L) / 1000000.0;

    if (print_) {
        outfile->Printf("    Estimated memory usage: %.1lf MB\n\n", memory);
    }
    if (options_.get_bool("SAPT_MEM_CHECK"))
        if (mem < vvnri + ovov * 3L) throw PsiException("Not enough memory", __FILE__, __LINE__);

    outfile->Printf("    Natural Orbital Cutoff: %11.3E\n", occ_cutoff_);
    // outfile->Printf("    Disp(T3) Truncation:    %11s\n", (nat_orbs_t3_ ? "Yes" : "No"));
    // outfile->Printf("    CCD (vv|vv) Truncation: %11s\n", (nat_orbs_v4_ ? "Yes" : "No"));
    outfile->Printf("    MBPT T2 Truncation:     %11s\n", (nat_orbs_t2_ ? "Yes" : "No"));
    outfile->Printf("\n");
}

void PSAPT::print_results() {
    if (e_ind30r_ != 0.0)
        e_exch_ind30r_ = e_ind30r_ * (e_exch_ind30_ / e_ind30_);
    else
        e_exch_ind30r_ = 0.0;
    
    double dHF2 = eHF_ - (e_elst10_ + e_exch10_ + e_ind20_ + e_exch_ind20_);
    // double dHF3 = eHF_ - (e_elst10_ + e_exch10_ + e_ind20_ +  e_exch_ind20_ + e_ind30r_ +
    //                            e_exch_ind30r_);

    e_sapt0_ = e_elst10_ + e_exch10_ + dHF2 + e_ind20_ + e_disp20_ + (e_exch_ind20_ + e_exch_disp20_);
    
    e_sapt2_ = e_elst10_ + e_exch10_ + dHF2 + e_ind20_ + e_disp20_ + (e_exch_ind20_ + e_exch_disp20_) +
               e_elst12_ + (e_exch11_ + e_exch12_ + e_exch_ind22_) + e_ind22_;

    double ind3 = e_ind30r_ + e_exch_ind30r_ ;
    double tot_elst = e_elst10_ + e_elst12_ + e_elst13_;
    double tot_exch = e_exch10_ + (e_exch11_ + e_exch12_);
    double tot_ind =  dHF2 + e_ind20_ + e_ind22_ + (e_exch_ind20_ + e_exch_ind22_ );
    double tot_ind3 = (ind3)  + e_ind20_ + e_ind22_ + (e_exch_ind20_ + e_exch_ind22_ );
    e_saptd4_ = tot_elst + tot_exch + tot_ind;

    double tot_ct = e_ind20_ + e_ind22_ + e_ind30r_ + (e_exch_ind20_ + e_exch_ind22_ );

    double tot_disp = 0.0;

    std::string scaled ="   ";
    outfile->Printf(
        "  "
        "--------------------------------------------------------------------------------------------------------"
        "\n");
    outfile->Printf("    Electrostatics            %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    tot_elst * 1000.0, tot_elst * pc_hartree2kcalmol, tot_elst * pc_hartree2kJmol);
    outfile->Printf("      Elst10,r                %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    e_elst10_ * 1000.0, e_elst10_ * pc_hartree2kcalmol, e_elst10_ * pc_hartree2kJmol);
    outfile->Printf("      Elst12,r                %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    e_elst12_ * 1000.0, e_elst12_ * pc_hartree2kcalmol, e_elst12_ * pc_hartree2kJmol);
    outfile->Printf("      Elst13,r                %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    e_elst13_ * 1000.0, e_elst13_ * pc_hartree2kcalmol, e_elst13_ * pc_hartree2kJmol);
    outfile->Printf("\n    Exchange %3s              %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), tot_exch * 1000.0, tot_exch * pc_hartree2kcalmol, tot_exch * pc_hartree2kJmol);
    outfile->Printf("      Exch10                  %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    e_exch10_ * 1000.0, e_exch10_ * pc_hartree2kcalmol, e_exch10_ * pc_hartree2kJmol);
    // outfile->Printf("      Exch10(S^2)             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
    //                 e_exch10_s2_ * 1000.0, e_exch10_s2_ * pc_hartree2kcalmol, e_exch10_s2_ * pc_hartree2kJmol);
    outfile->Printf("      Exch11(S^2) %3s         %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), e_exch11_ * 1000.0, e_exch11_ * pc_hartree2kcalmol,
                    e_exch11_ * pc_hartree2kJmol);
    outfile->Printf("      Exch12(S^2) %3s         %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), e_exch12_ * 1000.0, e_exch12_ * pc_hartree2kcalmol,
                    e_exch12_ * pc_hartree2kJmol);
    outfile->Printf("\n    Induction %3s             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), tot_ind * 1000.0, tot_ind * pc_hartree2kcalmol, tot_ind * pc_hartree2kJmol);
    outfile->Printf("\n    Induction(3) %3s             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), tot_ind3 * 1000.0, tot_ind3 * pc_hartree2kcalmol, tot_ind3 * pc_hartree2kJmol);
    outfile->Printf("      Ind20,r                 %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    e_ind20_ * 1000.0, e_ind20_ * pc_hartree2kcalmol, e_ind20_ * pc_hartree2kJmol);
    
    outfile->Printf("      Ind22                   %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    e_ind22_ * 1000.0, e_ind22_ * pc_hartree2kcalmol, e_ind22_ * pc_hartree2kJmol);
    outfile->Printf("      Exch-Ind20,r %3s        %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), e_exch_ind20_ * 1000.0,
                    e_exch_ind20_ * pc_hartree2kcalmol, e_exch_ind20_ * pc_hartree2kJmol);

    outfile->Printf("      Exch-Ind22 %3s          %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), e_exch_ind22_ * 1000.0,
                    e_exch_ind22_ * pc_hartree2kcalmol, e_exch_ind22_ * pc_hartree2kJmol);
    outfile->Printf("      delta HF,r (2) %3s      %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), dHF2 * 1000.0, dHF2 * pc_hartree2kcalmol, dHF2 * pc_hartree2kJmol);
    // outfile->Printf("      delta HF,r (3) %3s      %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
    //                 scaled.c_str(), dHF3 * 1000.0, dHF3 * pc_hartree2kcalmol, dHF3 * pc_hartree2kJmol);
    outfile->Printf("      Ind30,r                 %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                        e_ind30r_ * 1000.0, e_ind30r_ * pc_hartree2kcalmol, e_ind30r_ * pc_hartree2kJmol);

    outfile->Printf("      Exch-Ind30,r %3s        %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                            scaled.c_str(),  e_exch_ind30r_ * 1000.0,
                             e_exch_ind30r_ * pc_hartree2kcalmol,
                             e_exch_ind30r_ * pc_hartree2kJmol);
    outfile->Printf("\n    Dispersion %3s            %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), tot_disp * 1000.0, tot_disp * pc_hartree2kcalmol, tot_disp * pc_hartree2kJmol);
    

    outfile->Printf("\n  Total HF                    %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    eHF_ * 1000.0, eHF_ * pc_hartree2kcalmol, eHF_ * pc_hartree2kJmol);
    outfile->Printf("  Total SAPT0 %3s             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), e_sapt0_ * 1000.0, e_sapt0_ * pc_hartree2kcalmol, e_sapt0_ * pc_hartree2kJmol);
    outfile->Printf("  Total SAPT2 %3s             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scaled.c_str(), e_sapt2_ * 1000.0, e_sapt2_ * pc_hartree2kcalmol, e_sapt2_ * pc_hartree2kJmol);

    
    outfile->Printf("    delta(Ind - Ind(3))            %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    tot_ind - tot_ind3 * 1000.0, tot_ind - tot_ind3 * pc_hartree2kcalmol, tot_ind - tot_ind3 * pc_hartree2kJmol);    



        Process::environment.globals["SAPT-D4 ELST ENERGY"] = tot_elst;
        Process::environment.globals["SAPT-D4 IND ENERGY"] = tot_ind;
        Process::environment.globals["SAPT-D4 DISP ENERGY"] = 0.0;
        Process::environment.globals["SAPT-D4 EXCH ENERGY"] = tot_exch;
        Process::environment.globals["SAPT-D4 TOTAL ENERGY"] = tot_exch+tot_ind+tot_elst;

            // for auto-docs
            /*- Process::environment.globals["SAPT ELST ENERGY"] -*/
            /*- Process::environment.globals["SAPT EXCH ENERGY"] -*/
            /*- Process::environment.globals["SAPT IND ENERGY"] -*/
            /*- Process::environment.globals["SAPT DISP ENERGY"] -*/
            /*- Process::environment.globals["SAPT TOTAL ENERGY"] -*/

        
    // }
}
}  // namespace sapt
}  // namespace psi
