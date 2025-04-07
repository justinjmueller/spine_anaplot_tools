python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_1p_contained.root --config configurations/analyses/sbnd/1D/muon2024_1mu1p_contained_1D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_1p_allmuons.root --config configurations/analyses/sbnd/1D/muon2024_1mu1p_allmuons_1D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_Np_contained.root --config configurations/analyses/sbnd/1D/muon2024_1muNp_contained_1D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_Np_allmuons.root --config configurations/analyses/sbnd/1D/muon2024_1muNp_allmuons_1D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_X_contained.root --config configurations/analyses/sbnd/1D/muon2024_1muX_contained_1D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_X_allmuons.root --config configurations/analyses/sbnd/1D/muon2024_1muX_allmuons_1D.toml

python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_1p_contained.root --config configurations/analyses/sbnd/2D/muon2024_1mu1p_contained_2D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_1p_allmuons.root --config configurations/analyses/sbnd/2D/muon2024_1mu1p_allmuons_2D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_Np_contained.root --config configurations/analyses/sbnd/2D/muon2024_1muNp_contained_2D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_Np_allmuons.root --config configurations/analyses/sbnd/2D/muon2024_1muNp_allmuons_2D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_X_contained.root --config configurations/analyses/sbnd/2D/muon2024_1muX_contained_2D.toml
python3 spineplot.py --input ../cafana/output/larcv_sbnd_bnb_cosmics_spine_updated/larcv_sbnd_bnb_cosmics_spine_updated_X_allmuons.root --config configurations/analyses/sbnd/2D/muon2024_1muX_allmuons_2D.toml

mv ../output/1mu1p/* ../output/larcv_sbnd_bnb_cosmics_spine_updated/1mu1p/
mv ../output/1muNp/* ../output/larcv_sbnd_bnb_cosmics_spine_updated/1muNp/
mv ../output/1muX/* ../output/larcv_sbnd_bnb_cosmics_spine_updated/1muX/

