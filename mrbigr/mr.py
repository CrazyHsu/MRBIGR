#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import pyranges as pr
from sklearn.linear_model import LinearRegression
from pandas_plink import read_plink1_bin
from scipy.stats import chi2
import mrbigr.multiprocess as mu
import warnings
import shutil
import glob
import sys
import os


def lm_res(y, X):
    idx = np.isnan(X) | np.isnan(y)
    X = X[~idx].reshape(-1, 1)
    y = y[~idx]
    lm = LinearRegression().fit(X, y)
    sigma2 = np.sum((y - lm.predict(X))**2)/(X.shape[0]-1)
    se = np.sqrt(np.diag(np.linalg.pinv(np.dot(X.T, X))))[-1] * sigma2
    effect = lm.coef_[-1]
    rsq = lm.score(X, y)
    return pd.Series(dict(zip(['effect', 'se', 'rsq'], [effect, se, rsq])))


def MR(mTrait, pTrait, g, pvalue_cutoff):
    p = pd.concat([mTrait.to_frame(), pTrait], axis=1)
    snp_lm = p.apply(lm_res, args=[g.values])
    mTrait_lm = p.apply(lm_res, args=[mTrait.values])
    pTrait_mTrait = snp_lm.loc['effect', :][1:] / snp_lm.loc['effect', :][0]
    var_upper = pTrait.var() * (1 - mTrait_lm.loc['rsq', :][1:])
    var_down = mTrait.shape[0] * mTrait.var() * snp_lm.loc['rsq', :][0]
    var = var_upper / var_down
    TMR = pTrait_mTrait**2 / var
    pvalue = 1 - chi2.cdf(TMR, 1)
    MR_res = pd.DataFrame(dict(zip(['snp', 'mTrait', 'pTrait', 'effect', 'TMR', 'pvalue'],
                        [[g.name]*pTrait.shape[1], [mTrait.name]*pTrait.shape[1], pTrait.columns, pTrait_mTrait, TMR, pvalue])))
    MR_res = MR_res.loc[(MR_res.pvalue<=pvalue_cutoff) & (MR_res.mTrait!=MR_res.pTrait),:]
    return MR_res


def MR_parallel(mTrait_qtl, mTrait, pTrait, geno, threads, pvalue_cutoff):
    args = list()
    mTrait_qtl = mTrait_qtl.loc[mTrait_qtl.phe_name.isin(set(mTrait_qtl.phe_name) & set(mTrait.columns)), :]
    for index, row in mTrait_qtl.iterrows():
        rs = row['SNP']
        mTrait_name = row['phe_name']
        args.append((mTrait.loc[:, mTrait_name], pTrait, geno.loc[:, rs], pvalue_cutoff))
    res = mu.parallel(MR, args, threads)
    MR_res = pd.concat([i for i in res])
    return MR_res


def var_bXY(bzx, bzx_se, bzy, bzy_se):
    varbXY = bzx_se**2*bzy**2/bzx**4 + bzy_se**2/bzx**2
    return varbXY


def generate_geno_batch(mTrait_qtl, mTrait, pTrait, geno, threads, bed_dir, rs_dir):
    if os.path.exists(bed_dir):
        shutil.rmtree(bed_dir)
    os.mkdir(bed_dir)
    if os.path.exists(rs_dir):
        shutil.rmtree(rs_dir)
    os.mkdir(rs_dir)
    plink_extract = 'plink -bfile {} -extract {} --make-bed -out {}'
    geno_batch = list()
    mTrait_qtl = mTrait_qtl.loc[mTrait_qtl.phe_name.isin(set(mTrait_qtl.phe_name) & set(mTrait.columns)), :]
    for mTrait_name in mTrait_qtl.phe_name.unique():
        out_name = bed_dir.strip('/') + '/' + mTrait_name
        rs = mTrait_qtl.loc[mTrait_qtl.phe_name == mTrait_name, 'SNP']
        rs_name = rs_dir.strip('/') + '/' + '_'.join([mTrait_name,'rs.txt'])
        pd.Series(rs).to_frame().to_csv(rs_name, index=False, header=None)
        geno_batch.append((plink_extract.format(geno, rs_name, out_name),))
    out_name = bed_dir.strip('/') + '/pTrait'
    rs_name = rs_dir.strip('/') + '/pTrait_rs.txt'
    mTrait_qtl['SNP'].to_frame().to_csv(rs_name, index=False, header=None)
    geno_batch.append((plink_extract.format(geno, rs_name, out_name),))
    mu.parallel(mu.run, geno_batch, threads)
    for fn in glob.glob(bed_dir.strip('/')+'/*fam'):
        fam = pd.read_csv(fn, sep=' ', header=None)
        mTrait_name = fn.split('/')[-1].replace('.fam', '')
        if mTrait_name == 'pTrait':
            pTrait = pTrait.reindex(fam[0])
            fam.index = fam[0]
            fam = pd.concat([fam, pTrait], axis=1)
        else:
            fam.loc[:, 5] = mTrait.loc[:, mTrait_name].reindex(fam[0]).values
        fam.to_csv(fn, index=False, header=None, sep=' ', na_rep='NA')


def calc_MLM_effect(bed_dir, pTrait, threads, geno):
    args = list()
    geno_prefix = geno.split('/')[-1]
    fam = pd.read_csv(geno + '.fam', sep=r'\s+', header=None)
    fam[5] = 1
    fam.to_csv(geno_prefix + '.link.fam', sep='\t', na_rep='NA', header=None, index=False)
    if os.path.exists(geno_prefix + '.link.bed'):
        os.remove(geno_prefix + '.link.bed')
    if os.path.exists(geno_prefix + '.link.bim'):
        os.remove(geno_prefix + '.link.bim')
    os.symlink(geno + '.bed', geno_prefix + '.link.bed')
    os.symlink(geno + '.bim', geno_prefix + '.link.bim')
    related_matrix_cmd = 'gemma.linux -bfile {0}.link -gk 1 -o {1}'.format(geno_prefix, geno_prefix)
    s = mu.run(related_matrix_cmd)
    if s != 0:
        return None
    gemma_cmd_mTrait = 'gemma.linux -bfile {0} -k ./output/{1}.cXX.txt -lmm -n 1 -o {2}'
    gemma_cmd_pTrait = 'gemma.linux -bfile {0} -k ./output/{1}.cXX.txt -lmm -n {2} -o {3}'
    for i in glob.glob(bed_dir + '/*.bed'):
        i = i.replace('.bed', '')
        if i.split('/')[-1] != 'pTrait':
            prefix = i.split('/')[-1]
            args.append((gemma_cmd_mTrait.format(i, geno_prefix, 'mTrait_' + prefix),))
        else:
            for _, pTrait_name in enumerate(pTrait.columns):
                args.append((gemma_cmd_pTrait.format(i, geno_prefix, _ + 2, 'pTrait_' + pTrait_name),))
    s = mu.parallel(mu.run, args, threads)
    os.remove(geno_prefix + '.link.bed')
    os.remove(geno_prefix + '.link.bim')
    os.remove(geno_prefix + '.link.fam')
    return s


def get_MLM_effect(fn):
    assoc = pd.read_csv(fn, sep='\t')
    assoc.index = assoc['rs']
    return assoc[['beta', 'se']]


def get_MLM_effect_parallell(assoc_dir, mTrait, pTrait, threads):
    mTrait_effect = pd.DataFrame()
    args = []
    #pTrait_name = []
    # for fn in glob.glob(assoc_dir.strip('/') + '/mTrait*.assoc.txt'):
    #     mTrait_name = fn.split('/')[-1].split('_')[-1].replace('.assoc.txt', '')
    #     assoc = pd.read_csv(fn, sep='\t')
    #     assoc.index = mTrait_name+';' + assoc['rs']
    #     mTrait_effect = pd.concat([mTrait_effect, assoc[['beta', 'se']]])
    # for fn in glob.glob(assoc_dir.strip('/') + '/pTrait*assoc.txt'):
    #     pTrait_name.append(fn.split('/')[-1].split('_')[-1].replace('.assoc.txt', ''))
    #     args.append((fn,))
    for mTrait_name in mTrait.columns:
        fn = assoc_dir.strip('/') + '/mTrait_' + mTrait_name + '.assoc.txt'
        assoc = pd.read_csv(fn, sep='\t')
        assoc.index = mTrait_name + ';' + assoc['rs']
        mTrait_effect = pd.concat([mTrait_effect, assoc[['beta', 'se']]])
    for pTrait_name in pTrait.columns:
        fn = assoc_dir.strip('/') + '/pTrait_' + pTrait_name + '.assoc.txt'
        args.append((fn,))
    pTrait_res = mu.parallel(get_MLM_effect, args, threads)
    pTrait_effect = pd.concat([i['beta'] for i in pTrait_res], axis=1)
    pTrait_effect.columns = pTrait.columns
    pTrait_se = pd.concat([i['se'] for i in pTrait_res], axis=1)
    pTrait_se.columns = pTrait.columns
    return mTrait_effect, pTrait_effect, pTrait_se


def MR_MLM(mTrait_effect_snp, pTrait_effect_snp, pTrait_se_snp, pvalue_cutoff):
    mTrait_name, rs = mTrait_effect_snp.name.split(';')
    bxy = pTrait_effect_snp / mTrait_effect_snp['beta']
    varbXY = var_bXY(mTrait_effect_snp['beta'], mTrait_effect_snp['se'], pTrait_effect_snp, pTrait_se_snp)
    TMR = bxy**2 / varbXY
    pvalue = 1 - chi2.cdf(TMR, 1)
    MR_res = pd.DataFrame(dict(zip(['snp', 'mTrait', 'pTrait', 'effect', 'TMR', 'pvalue'],
                                   [[rs] * pTrait_effect_snp.shape[0], [mTrait_name] * pTrait_effect_snp.shape[0],
                                    pTrait_effect_snp.index, bxy, TMR, pvalue])))
    MR_res = MR_res.loc[(MR_res.pvalue<=pvalue_cutoff) & (MR_res.mTrait!=MR_res.pTrait), :]
    return MR_res


def MR_MLM_parallel(mTrait_qtl, mTrait_effect, pTrait_effect, pTrait_se, threads, pvalue_cutoff):
    args = []
    mTrait_qtl = mTrait_qtl.loc[mTrait_qtl.phe_name.isin(set(mTrait_qtl.phe_name) & set(mTrait_effect.index.map(lambda x:x.split(';')[0]))), :]
    for index, row in mTrait_qtl.iterrows():
        mTrait_name = row['phe_name']
        rs = row['SNP']
        args.append((mTrait_effect.loc[';'.join([mTrait_name, rs]),:], pTrait_effect.loc[rs, :], pTrait_se.loc[rs, :], pvalue_cutoff))
    res = mu.parallel(MR_MLM, args, threads)
    MR_res = pd.concat([i for i in res])
    return MR_res


def target_type(qtl, tf, target):
    tf_qtl = qtl.loc[qtl.phe_name.isin(tf.gene_id), :]
    target_qtl = qtl.loc[qtl.phe_name.isin(target.gene_id), :]
    target_qtl_peak_pos = target_qtl[['CHR', 'phe_name']]
    target_qtl_peak_pos.loc[:, 'start'] = target_qtl['SNP'].apply(lambda x: int(x.split('_')[-1]) - 1)
    target_qtl_peak_pos.loc[:, 'end'] = target_qtl['SNP'].apply(lambda x: int(x.split('_')[-1]) + 1)
    target_qtl_peak_pos = target_qtl_peak_pos[['CHR', 'start', 'end', 'phe_name']]
    target_qtl_peak_pos.columns = ['Chromosome', 'Start', 'End', 'phe_name']
    tf_qtl_peak_range = tf_qtl[['CHR', 'phe_name']]
    tf_qtl_peak_range.loc[:, 'start'] = tf_qtl['SNP'].apply(lambda x: int(x.split('_')[-1]) - 500000)
    tf_qtl_peak_range.loc[:, 'end'] = tf_qtl['SNP'].apply(lambda x: int(x.split('_')[-1]) + 500000)
    tf_qtl_peak_range = tf_qtl_peak_range[['CHR', 'start', 'end', 'phe_name']]
    tf_qtl_peak_range.columns = ['Chromosome', 'Start', 'End', 'phe_name']
    target_qtl_peak_pos = pr.PyRanges(target_qtl_peak_pos)
    tf_qtl_peak_range = pr.PyRanges(tf_qtl_peak_range)
    coloc = tf_qtl_peak_range.join(target_qtl_peak_pos)
    coloc_df = pd.DataFrame()
    for k in sorted(coloc.dfs.keys()):
        coloc_df = pd.concat([coloc_df, coloc.dfs[k]])
    return coloc_df


def read_genotype(geno_prefix):
    try:
        G = read_plink1_bin(geno_prefix + '.bed', geno_prefix + '.bim', geno_prefix + '.fam', ref='a0',verbose=False)
    except Exception:
        return None
    return G
