c
c   $Id: olemcdf.inc,v 1.2 1999/08/16 20:22:32 eds Exp $
c
c common for cdf file for fields
c
      integer fldid(5),ufid(5),vfid(5),wfid(5),pfid(5),
     *     thfid(5),salfid(5),tkefid(5),rho_pfid(5),khfid(5),
     *     kmfid(5),rhofid(5),hf_sgsid(5),strtf(4),countf(4),
     *     tindexf(5),countuf(4),countvf(4),countwf(4),tfid(5),
     *     trcrfid(5,trcdim),zifid(5),w_starfid(5),t_starfid(5),
     *     qfid(5),u_starfid(5),hf_topfid(5),ustr_tfid(5),
     *     vstr_tfid(5),wave_lfid(5),wave_hfid(5),w_anglefid(5),
     *     swf_topfid(5),rainfid(5),lhf_topfid(5)
c
      common /cdffld/fldid,ufid,vfid,wfid,pfid,thfid,salfid,tkefid,
     *     rho_pfid,trcrfid,khfid,kmfid,rhofid,
     *     hf_sgsid,strtf,countf,countuf,countvf,countwf,tfid,tindexf,
     *     zifid,w_starfid,t_starfid,qfid,u_starfid,
     *     hf_topfid,ustr_tfid,vstr_tfid,wave_lfid,wave_hfid,
     *     w_anglefid,swf_topfid,rainfid,lhf_topfid
c
c common for cdf file for time series
c
      integer timid,utid,vtid,wtid,ptid,thtid,saltid,tketid,
     *     rho_ptid,trcrtid(trcdim),khtid,kmtid,ttid,tindext,
     *     t_bartid,s_bartid,zitid,w_startid,t_startid,qtid,
     *     u_startid,t_avetid,t2_avetid,s_avetid,s2_avetid,
     *     tke_avetid,KE_avetid,w2_avetid,wt_avetid,ws_avetid,
     *     wrho_avetid,rhop_avetid,u_avetid,v_avetid,u2_avetid,
     *     v2_avetid,hf_avetid,rhof_avetid,sf_avetid,p2_avetid,
     *     ufilt_avetid,vfilt_avetid,wfilt_avetid,p_avetid,
     *     a_avetid,b_avetid,sd_avetid,sg_avetid,uw_avetid,
     *     vw_avetid,sp_avetid,hf_toptid,ustr_ttid,vstr_ttid,
     *     wave_ltid,wave_htid,w_angletid,swf_toptid,raintid,
     *     lhf_toptid,S_0tid,L_atid,km_avetid,petid,pettid,pewtid,
     *     sig_avetid,disptid,cfratid,wctid,hbctid,chtid,sftid,sbctid,
     *     wtatid,wsatid,bvsttid,wvsttid,uvdutid,uvdvtid,uwdutid,
     *     vwdvtid,wrhotid
c
      common /cdftim/timid,utid,vtid,wtid,ptid,thtid,saltid,tketid,
     *     rho_ptid,trcrtid,khtid,kmtid,ttid,tindext,
     *     t_bartid,s_bartid,zitid,w_startid,t_startid,qtid,
     *     u_startid,t_avetid,t2_avetid,s_avetid,s2_avetid,
     *     tke_avetid,KE_avetid,w2_avetid,wt_avetid,ws_avetid,
     *     wrho_avetid,rhop_avetid,u_avetid,v_avetid,u2_avetid,
     *     v2_avetid,hf_avetid,rhof_avetid,sf_avetid,p2_avetid,
     *     ufilt_avetid,vfilt_avetid,wfilt_avetid,p_avetid,
     *     a_avetid,b_avetid,sd_avetid,sg_avetid,uw_avetid,
     *     vw_avetid,sp_avetid,hf_toptid,ustr_ttid,vstr_ttid,
     *     wave_ltid,wave_htid,w_angletid,swf_toptid,raintid,
     *     lhf_toptid,S_0tid,L_atid,km_avetid,petid,pettid,pewtid,
     *     sig_avetid,disptid,cfratid,wctid,hbctid,chtid,sftid,sbctid,
     *     wtatid,wsatid,bvsttid,wvsttid,uvdutid,uvdvtid,uwdutid,
     *     vwdvtid,wrhotid
c
c common for cdf file of input time series
c
      integer bcid,hftid,ustrtid,vstrtid,wavelid,wavehid,angleid,bctid,
     *        swhftid,rainid,lhftid
c
      common /cdfbc/bcid,hftid,ustrtid,vstrtid,wavelid,wavehid,angleid,
     *       bctid,swhftid,rainid,lhftid
c
c common for cdf file for slabs
c
      integer ufidx,vfidx,wfidx,pfidx,thfidx,salfidx,tkefidx,
     *       ufidy,vfidy,wfidy,pfidy,thfidy,salfidy,tkefidy,rho_pfidx,
     *       ufidz,vfidz,wfidz,pfidz,thfidz,salfidz,tkefidz,hifidz,
     *       rho_pfidz,rho_pfidy,khfidx,khfidy,khfidz,kmfidx,kmfidy,
     *       kmfidz,rhofidx,rhofidy,rhofidz,hf_sgsfidx,hf_sgsfidy,
     *       hf_sgsfidz,trcrfidx,trcrfidy,trcrfidz,botfidz,hfifidz,
     *       uffidz,vffidz,wffidz,rffidz,vafidx,vafidy,vafidz,
     *       dcfidx,dcfidy,dcfidz,dpfidx,dpfidy,dpfidz,
     *       sgfidx,sgfidy,sgfidz,sfifidz,slbid,tsid,tindexs
c
      common /cdfslb/ufidx(5),vfidx(5),wfidx(5),pfidx(5),thfidx(5),
     *     salfidx(5),tkefidx(5),ufidy(5),vfidy(5),wfidy(5),pfidy(5),
     *     thfidy(5),salfidy(5),tkefidy(5),rho_pfidx(5),ufidz(5),
     *     vfidz(5),wfidz(5),pfidz(5),thfidz(5),salfidz(5),tkefidz(5),
     *     rho_pfidz(5),rho_pfidy(5),khfidx(5),khfidy(5),khfidz(5),
     *     kmfidx(5),kmfidy(5),kmfidz(5),rhofidx(5),rhofidy(5),
     *     rhofidz(5),hf_sgsfidx(5),hf_sgsfidy(5),hf_sgsfidz(5),
     *     trcrfidx(5,trcdim),trcrfidy(5,trcdim),trcrfidz(5,trcdim),
     *     hifidz(5),botfidz(5),slbid(5),tsid(5),tindexs(5),
     *     uffidz(5),vffidz(5),wffidz(5),rffidz(5),
     *     dcfidx(5),dcfidy(5),dcfidz(5),vafidx(5),vafidy(5),vafidz(5),
     *     dpfidx(5),dpfidy(5),dpfidz(5),sfifidz(5),hfifidz(5),
     *     sgfidx(5),sgfidy(5),sgfidz(5)

