# LayeredSpheres_PKPD_ranibizumab

## Overview
The combined mathematical model of drug release dynamics from bi-layered core-shell polymeric microspheres with the three-compartment pharmacokinetic/pharmacodynamic model of ranibizumab and VEGF binding to evaluate the biodistribution of ranibizumab and quantify its effects on VEGF inhibition in three compartments: retina, vitreous humor, and aqueous humor. The combined model was used to determine how initial drug loading and design parameters influence VEGF pharmacodynamic suppression duration for single- and bi-layered microspheres. 
## Authors
Md Tanben Rahman<sup>a</sup>, Mohammad Aminul Islam<sup>a</sup>, Koki Kanehira<sup>a</sup>, Sarita Das<sup>a</sup>, Yaman Oklla<sup>a</sup>,
Eduardo A. Chacin Ruiz<sup>a</sup>, Katelyn E. Swindle-Reilly<sup>e,f,g</sup>, Ashlee N. Ford Versypt<sup>a,b,c,d</sup><br/>

<sup>a</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY, 14260, USA<br/>
<sup>b</sup>Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, 14260, USA<br/>
<sup>c</sup>Institute for Artificial Intelligence and Data Science, University at Buffalo, The State University of New York, Buffalo, NY, 14260, USA<br/>
<sup>d</sup>Department of Pharmaceutical Sciences, University at Buffalo, The State University of New York, Buffalo, NY, 14215, USA<br/>
<sup>e</sup>William G. Lowrie Department of Chemical and Biomolecular Engineering, The Ohio State University, Columbus, OH, 43210, USA<br/>
<sup>f</sup>Department of Ophthalmology and Visual Sciences, The Ohio State University, Columbus, OH, 43210, USA<br/>
<sup>g</sup>Department of Biomedical Engineering, The Ohio State University, Columbus, OH, 43210, USA<br/>

## Manuscript

## Scripts

* **without_DDS.m** This file simulates the dynamics of ranibizumab and free VEGF and the times for the VEGF concentration to return to 10% and 50% of baseline concentration for different doses of ranibizumab without DDS

* **DDS_doses.m** This file simulates the dynamics of ranibizumab and free VEGF, the times for the VEGF concentration to return to 10% and 50% of baseline concentration, and drug release dynamics of bi-layered DDS for different doses of ranibizumab without DDS

* **bi_layer_chitosan.m** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration and drug release dynamics for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor. The radius of the chitosan layer is varied in the simulation

* **bi_layered_PCL.m** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration and drug release dynamics for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor. The radius of the PCL layer is varied in the simulation

* **bi_layer_changing_both.m** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration and drug release dynamics for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor. The radius of both chitosan and PCL layer is varied simultaneously in the simulation

* **DDS_k_D_k_on_variations.m** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor for varying dissociation constants $K_{D}$ with $k_{\text{off}}$ fixed at 0.864 day$^{-1}$

* **DDS_k_D_k_on_variations.m** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor for varying dissociation constants $K_{D}$ with $k_{\text{off}}$ fixed at 0.864 day$^{-1}$

* **without_DDS_k_D_k_on_variations.m** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg without DDS at retina, vitreous, and aqueous humor for varying dissociation constants $K_{D}$ with $k_{\text{off}}$ fixed at 0.864 day$^{-1}$

* **DDS_k_off_variations.m** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor for varying $k_{\text{off}}$ with constant $K_{D}$

* **without_DDS_k_off_variations.m** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg without DDS at retina, vitreous, and aqueous humor for varying $k_{\text{off}}$ with constant $K_{D}$

* **FD_spheres_variable_diffusivity_two_spheres.m** This is the function file to solve the PDE of bi-layered chitosan-PCL core-shell DDS for Fickian diffusion within a radially symmetric sphere

* **solve_FD_spheres_variable_diffusivity.m** This file solves the PDE of bi-layered chitosan-PCL core-shell DDS for Fickian diffusion within a radially symmetric sphere

## Acknowledgements
This work was supported by National Institutes of Health grant R35GM133763 to ANFV, R01EB032870 to KESR and ANFV, the Owen Locke Foundation to KESR, and the University at Buffalo. We thank lab members and committee member Dr. Rudiyanto Gunawan for their thorough feedback on this manuscript and helpful discussions.