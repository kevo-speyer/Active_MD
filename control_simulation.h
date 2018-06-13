/* Incorporates C preprocessor as a way of compiling the program with the appropriate physics */

/* SYSTEM: determines the system type compiled in the program
*  0= brush-melt channel, no PBC in z
 * 1= droplet(no brush on top wall) 
 * 2= charged system.No brush, 
 * 3= charged channel+brush  */

#define SYSTEM 1 


/*BRUSH TYPE: sets how the grafted beads from the polymer brush will be chosen
 * 0: randomly distributed with a uniform distribution. Beads may overlap one another.
 * 1: randomly distributed with a uniform distribution. Overlaping beads are forbidden.
 * 2: ordered brush. The distance between neighbouring grafted beads is equal for all points
 * 3: ordered brush, starts aligned in one direction, i.e. 'combed'*/

#define BRUSH_TYPE 3 

/* SYMMETRY */
/* 0 = channel geometry: no PBC in Z; 
 * 1 = Bulk geometry PBC in 3D */

#define SYMMETRY 0   

/* Flat dynamics. Sets every 'y' component to zero.*/

#undef BIDIMENSIONAL

/* Potentials and particles */

/* Wall type 
 * 1: explicit wall 
 * 2: implicit wall 9-3, 
 * 3: top: hard wall, bottom: 9-3 4=top and bottom, hard walls   (droplets)
 * 4: hard walls in top and bottom walls */

#define WALL 3

/*      Potential cut-offs for LJ
 * 3= non-additive+poor 
 * 2= non-additive+good 
 * 1= good solvent 
 * 0= poor solvent      */

/* NOTE:  WALL is  not used if SYMMETRY /= 1  */
#define SOLVENT 1    

#undef HYDROPHOBIA  /* if def, the interaction between brush and melt is purely repulsive *//
#undef BRUSH_IN_GOOD_SV /*if def the interaction between grafted polimers is purely repulsive*/
#define BENDING  /* if def the grafted polymers are assumed to be semiflxible: bending potential *//
#undef BENDING_MELT  /* if def the melt polymers are assumed to be semiflxible: bending potential *//
#define ORIENTATION  /* if def the grafted polymers will be oriented through an harmonic potential*//

#undef PARTICLE_4 /* If defined the program runs with four different particle type */

#undef STARS /*whether you want to simulate with or without stars, sigma is fixed to 1. As well as sigma for walls */

/* Active brush polymers*/

#define ACTIVE_BRUSH

#define METRO 1 /* 0 is 2D activation model. Angle meassured in x dimension*/
                /* 1 is 3D act. model. Activation angle is meassured in the plane (isotropic) */
/* Brushes are coupled by springs */

#undef SPRING_ARRAY

#undef RANGLE

/* Thermostat */

#define THERMOSTAT 0 /*  1=LGV 0=DPD       */

#define DPD_WEIGHT 0 /*  0=usual choice of DPD weight: Wd=(1-r/rc)^2 ; 1= constant: Wr=Wd=1 ; 2 "quartic" */ 
                     /*  wd=(1-r/rc)^4                                                                    */  
#define DPD_CUT_OFF  2.24 /* if defined, takes this cutoff for DPD forces */

#undef DPD_VV       /* Adds a re-calculation of the Fd at the end of the iteration cycle. 
                     * Improves T=const ? Vattullainen JCP 116 2967 
                     * Good for soft potentials (not implemmented) */ 

/*
* Relaxation mechanisms when thermalizing new configurations 
*/

#define FORCE_SWITCH_ON /* switchs on the real force progressively to avoid overflows when there are strong
                       * particles overlaps */
                    /*NOTE: Use with caution because the program will not be doing real MD */
#undef RELAX       /* When defined, velocities  are set to zero in each time step. 
                    * This is not MD, but force relaxation */


#undef POISEUILLE    /* Adds external constant force to simulate Poiseuille flow      */
#undef SHEARED      /* if defined, the shear protocols are applied, mfa_input is different!! */
                    /* NOTE: if it is not defined, wall velocities can anyway been used */

#define  STORE 0 /* 0=Writes out folded coordinates,1=writes out unfolded (to follow diffusion)   */
#define  WR_FORCES 0 /* 0=Writes out sum of forces on brush heads,1=writes out individual forcen on brush heads  */
#define  PROFILES      /* this uses mide for calculating profiles in the MD run.       */
#undef DIFF /*ifdef, the program will calculate diffusion coefficients. Not implemmented for all the systems  */
#undef BRUSH_SEP /* ifdef separates the density profiles of top and bottom brushes */                    
#undef FLUKT /*write out observations of the forces between Particle 4 and 3 in file  fort.555*/
#undef PINNED 2 /*applicable only with particle 4 and particle 3 and brushes,*/
                 /* 1 - fixed position, thermostat is not applied */
                 /* 2 - with the spring, thermostat is applied*/
                 /* the distance between them - from system input*/
/* Misc variables */
#define NO_WARNS     /* Does not print the warning messages of the "beads too close"                      */     
/*#define FLUID_ROUTINE 0*/  /* controls if it uses the normal fluid routine (0) or the experimental HPC-tuned (1) */
#undef FREE_HEADS        /* if defined, the heads of brushes are not fixed, epsilon between head and surface is 250, potential - attractive*/
#define BIN_TYPE 2 /*0: uses binning.f90; single counting of each interaction;*/
                    /*1: uses my_binning.f90, from S. Plimpton and Cem Servantie versions; double counts each interaction.*/
                    /*2: uses cell_list.f90; does cell-linked lists, No Verlet-List. by Kevo*/

#define L_BOX 1 /*Only relevant if bin_type=2, */
                /* 1 is length of binning box = r_cut_max. 13 neighbor cells*/
                /* 2 is length of binning box = r_cut_max / 2 . 62 neighbor cells */
 
#define RESPA /*Implement multiple time scale molecular dynamics (Ref: Tuckerman Berne Martyna, 1992, JCP)*/

#define SOL_SOL_INT 2 /* Control the interactions between solvent-solvent. WARNING: IT ONLY WORKS WITH RESPA DEFINED(a_type=3) */
                      /* 1 is Lennard-Jones Potential*/
                      /* 2 is Soft Potentail (typically used in DPD)*/

#define BRUSH_SOL_INT 2 /* Switch between Lennard-Jones and Soft Potential for brush monomer and solvent interaction */                
                        /* 1  is Lennard-Jones Potential*/
                        /* 2 is Soft Potentail (typically used in DPD)*/

#define VEL_INIT 1  /* 1 Set initial velocities according to Maxwell-Boltzmann distro*/
                    /* 0 Set initial velocities to 0 */
