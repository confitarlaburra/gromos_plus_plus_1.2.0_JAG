# This is the programs automake template.
AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src

EXTRA_DIST =   0index_contrib.doxy

bin_PROGRAMS = \
	aggregates\
	rottemp\
	ener_rplc\
	espmap\
	lb_top\
	m_ener\
	prep_bb\
	cg_potential\
	countbb\
	stacking\
	shake\
	model\
	gca_ener\
	renumber\
	distmat\
	iwdcf\
	bb2tex\
	crossdih_top\
	saltbridge\
	bilayer_gel\
	make_rdc_spec\
	svd_fit\
	add_atom\
	top_to_build\
	helix_axis\
	rdc_sh\
	gromacs2gromos\
	diff_top\
	matrix_overlap\
	grid_dens\
	pdb2seq\
	nscatt\
	structal\
	cgLJpot\
	hvap\
	depsi2native\
	rep_rewrite\
        local_comp\
        prep_hybrid\
	sasa_new\
	ener_test\
	force\
        mut_top\
        cavity_load\
        temperature\
	gatherNT\
	pore_loading\
	gatherNT\
	pore_ads\
	pore_ax_rad\
	pore_pf\
	pore_pf_multi\
	pore_diffus\
	pore_P1\
	pore_Pcorr\
	pore_rates\
	pore_dist_corr\
	pore_surv_P\
	pore_hsolv\
	pore_hbond\
	graphite\
	pore_rot_rel\
	dipole_corr\
	Kinetic_MC\
	diffus_test
aggregates_SOURCES = aggregates.cc              
rottemp_SOURCES = rottemp.cc
ener_rplc_SOURCES = ener_rplc.cc
espmap_SOURCES = espmap.cc
lb_top_SOURCES = lb_top.cc
m_ener_SOURCES = m_ener.cc
prep_bb_SOURCES = prep_bb.cc
cg_potential_SOURCES = cg_potential.cc
countbb_SOURCES = countbb.cc
stacking_SOURCES = stacking.cc
shake_SOURCES = shake.cc
model_SOURCES = model.cc
gca_ener_SOURCES = gca_ener.cc
renumber_SOURCES = renumber.cc
distmat_SOURCES = distmat.cc
iwdcf_SOURCES = iwdcf.cc
bb2tex_SOURCES = bb2tex.cc
crossdih_top_SOURCES = crossdih_top.cc
saltbridge_SOURCES = saltbridge.cc
bilayer_gel_SOURCES = bilayer_gel.cc
make_rdc_spec_SOURCES = make_rdc_spec.cc
svd_fit_SOURCES = svd_fit.cc
add_atom_SOURCES = add_atom.cc
top_to_build_SOURCES = top_to_build.cc
helix_axis_SOURCES = helix_axis.cc
rdc_sh_SOURCES = rdc_sh.cc
gromacs2gromos_SOURCES = gromacs2gromos.cc
diff_top_SOURCES = diff_top.cc
matrix_overlap_SOURCES = matrix_overlap.cc
grid_dens_SOURCES = grid_dens.cc
pdb2seq_SOURCES = pdb2seq.cc
nscatt_SOURCES = nscatt.cc
structal_SOURCES = structal.cc
cgLJpot_SOURCES = cgLJpot.cc
hvap_SOURCES = hvap.cc
depsi2native_SOURCES = depsi2native.cc
rep_rewrite_SOURCES = rep_rewrite.cc
local_comp_SOURCES = local_comp.cc
prep_hybrid_SOURCES = prep_hybrid.cc
sasa_new_SOURCES = sasa_new.cc
ener_test_SOURCES = ener_test.cc
force_SOURCES = force.cc
mut_top_SOURCES = mut_top.cc
cavity_load_SOURCES = cavity_load.cc
temperature_SOURCES = temperature.cc
gatherNT_SOURCES =    gatherNT.cc
pore_loading_SOURCES = pore_loading.cc
gatherNT_SOURCES = gatherNT.cc
pore_ads_SOURCES = pore_ads.cc
pore_ax_rad_SOURCES = pore_ax_rad.cc
pore_pf_SOURCES = pore_pf.cc
pore_pf_multi_SOURCES = pore_pf_multi.cc
pore_diffus_SOURCES= pore_diffus.cc
pore_P1_SOURCES= pore_P1.cc
pore_Pcorr_SOURCES=pore_Pcorr.cc
pore_rates_SOURCES=pore_rates.cc
pore_dist_corr_SOURCES=pore_dist_corr.cc
pore_surv_P_SOURCES=pore_surv_P.cc
pore_hsolv_SOURCES=pore_hsolv.cc
pore_hbond_SOURCES=pore_hbond.cc
graphite_SOURCES=graphite.cc
pore_rot_rel_SOURCES=pore_rot_rel.cc
dipole_corr_SOURCES=dipole_corr.cc
Kinetic_MC_SOURCES=Kinetic_MC.cc
diffus_test_SOURCES=diffus_test.cc
LDADD = $(top_builddir)/src/libgromos.la

shake_LDADD = $(MDPP_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(MDPP_LDFLAGS)

gromacs2gromos_LDADD = $(GMX_LIB) $(top_builddir)/src/libgromos.la
gromacs2gromos_LDFLAGS = $(GMX_LDFLAGS) $(top_builddir)/src/libgromos.la
