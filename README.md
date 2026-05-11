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

* **Ranibizumab_DDS_All_Case.m** This file simulates the dynamics of ranibizumab and free VEGF concentrations in the retina, vitreous, and aqueous compartments. We evaluated the impact on the time required to return to 10% and 50% of baseline levels, considering the interplay between dose, microsphere type and size, and binding kinetics under different case scenarios. The following case scenarios were considered for evaluation

  **without_DDS_dose** This case simulates the dynamics of ranibizumab and free VEGF and the times for the VEGF concentration to return to 10% and 50% of baseline concentration for different doses of ranibizumab without DDS

  **With_DDS_dose** This case simulates the dynamics of ranibizumab and free VEGF, the times for the VEGF concentration to return to 10% and 50% of baseline concentration, and drug release dynamics of bi-layered DDS for different doses of ranibizumab with DDS

  **Chitosan_single** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration and drug release dynamics for a dose of 2 mg from single-layered chitosan core DDS at retina, vitreous, and aqueous humor. The radius of the chitosan is varied in the simulation

  **PCL_single** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration and drug release dynamics for a dose of 2 mg from single-layered PCL core DDS at retina, vitreous, and aqueous humor. The radius of the PCL is varied in the simulation

  **BiLayer_Changing_Chitosan** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration and drug release dynamics for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor. The radius of the chitosan layer is varied in the simulation

  **BiLayer_Changing_PCL** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration and drug release dynamics for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor. The radius of the PCL layer is varied in the simulation

  **BiLayer_changing_Both** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration and drug release dynamics for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor. The radius of both chitosan and PCL layer is varied simultaneously in the simulation

  **Without_DDS_changing_kD** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg without DDS at retina, vitreous, and aqueous humor for varying dissociation constants $K_{D}$ with $k_{\text{off}}$ fixed at 0.864 day$^{-1}$

  **With_DDS_changing_kD** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor for varying dissociation constants $K_{D}$ with $k_{\text{off}}$ fixed at 0.864 day$^{-1}$

  **DDS_k_D_k_on_variations** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor for varying dissociation constants $K_{D}$ with $k_{\text{off}}$ fixed at 0.864 day$^{-1}$

  **Without_DDS_changing_koff** This case simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg without DDS at retina, vitreous, and aqueous humor for varying $k_{\text{off}}$ with constant $K_{D}$

  **With_DDS_changing_Koff** This file simulates the times for the VEGF concentration to return to 10% and 50% of baseline concentration for a dose of 2 mg from bi-layered chitosan-PCL core-shell DDS at retina, vitreous, and aqueous humor for varying $k_{\text{off}}$ with constant $K_{D}$

* **FD_spheres_variable_diffusivity_two_spheres.m** This is the function file to solve the PDE of bi-layered chitosan-PCL core-shell DDS for Fickian diffusion within a radially symmetric sphere

* **solve_FD_spheres_variable_diffusivity.m** This file solves the PDE of bi-layered chitosan-PCL core-shell DDS for Fickian diffusion within a radially symmetric sphere
  
* **ODEs.m** This function file defines the system of ordinary differential equations used to simulate the dynamics of ranibizumab, free VEGF, VEGF–ranibizumab complex, and ranibizumab–VEGF–ranibizumab complex in the retina, vitreous, and aqueous humor compartments. It incorporates binding and unbinding kinetics, intercompartmental transport, elimination from the vitreous and aqueous humor, endogenous VEGF production, and drug input from the DDS release profile
  
* **ScriptForExportingImages.m** exports standardized images in .pdf and .tiff format
  
* **All the .mat** files contain the information generated from the different case scenarios that are used for creating the figures and Table
  
*	**Table_Figure_together_Ranibizumab.m** This file generates a .csv file that combines all data and includes all necessary figures
  
*	**Patient40Data.xlsx** This file contains the patient data on free VEGF concentration



## Acknowledgements
This work was supported by National Institutes of Health grant R35GM133763 to ANFV, R01EB032870 to KESR and ANFV, the Owen Locke Foundation to KESR, and the University at Buffalo. We thank lab members and committee member Dr. Rudiyanto Gunawan for their thorough feedback on this manuscript and helpful discussions.
