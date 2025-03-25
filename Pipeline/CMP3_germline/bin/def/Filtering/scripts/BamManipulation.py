#! /usr/bin/env python3
import pysam
import numpy as np
import matplotlib.pyplot as mp
import matplotlib.gridspec as gridspec
"""
>>> import pysam
>>> samfile=pysam.AlignmentFile("/home/flanduzzi@iit.local/Downloads/mapt.NA12156.altex.bam")
>>> print(samfile.header)
>>> print(samfile.check_index(())
>>> samfile.get_index_statistics()
[IndexStats(contig='NC_000017.10', mapped=0, unmapped=0, total=0), IndexStats(contig='NT_167251.1', mapped=0, unmapped=0, total=0)]
>>> for rr in samfile.fetch('NC_000017.10', 20, 30) :
...     print(rr)
"""

def BINS_rice_rule(n):
    return int(np.ceil(2*np.power(n,(1./3.))))

def BINS_struges_rule(n):
    return int(np.log2(n)+1.)

def find_binsize(num_points, target_coverage, binning_rules=BINS_rice_rule):
    return int(max(np.ceil(target_coverage/2./BINS_rice_rule(num_points)),1))

class BAM_Manager():

    def __init__(self, source_file):
        self.bam=pysam.AlignmentFile(source_file, mode="rb")
        self.verbose=False

    def histogram(self, contig_name, start, stop, bin_size=1, min_coverage=0, max_coverage=5000, normalized=False):
        hist_elements=int(np.ceil((max_coverage-min_coverage+0.0)/bin_size))
        hist_cov=np.zeros(hist_elements,dtype=int)
        average_in_range=0.
        sd_in_range=0.
        maxim=stop-start
        bases_ouside_range=0
        ### CHECK each position on the BAM
        for el in self.bam.pileup(contig_name, start, stop):
            pos=el.pos-start
            #### CHECK that POSITION is really inside the rage
            if pos>=0 and pos <= maxim :
                num_read=el.n
                ### CHECK if the coverage is inside the BOUNDARIES
                if num_read>=min_coverage and num_read<max_coverage:
                    average_in_range+=num_read
                    sd_in_range+=num_read*num_read
                    hist_cov[int((num_read-min_coverage)/bin_size)]+=1
            else :
                bases_ouside_range+=1

        print("Bases outside RANGE:", bases_ouside_range)
        bases_in_range=np.sum(hist_cov)
        if normalized and bases_in_range>0 :
            hist_cov=hist_cov/(bases_in_range+0.0)
        sd_in_range=np.sqrt((sd_in_range - average_in_range*average_in_range/bases_in_range)/bases_in_range) if bases_in_range>0 else 0.0
        if bases_in_range>0 :
            average_in_range/=bases_in_range
        if self.verbose:
            print(f"Region {contig_name},{start}:{stop} [{min_coverage}:{max_coverage}] bases:{bases_in_range} avg:{average_in_range} stddev:{sd_in_range}")
        return hist_cov, bases_in_range, average_in_range, sd_in_range

    def histogram_and_coverage(self, contig_name, start, stop, bin_size=1, min_coverage=0, max_coverage=5000, normalized=False):
        hist_elements=int(np.ceil((max_coverage-min_coverage+0.0)/bin_size))
        hist_cov=np.zeros(hist_elements,dtype=int)
        average_in_range=0.
        sd_in_range=0.

        maxim=stop-start
        bases_ouside_range=0
        ### CREATE the EMPTY vector for the CURVE
        data=np.zeros((2,maxim+1), dtype=float)
        ### CHECK each position on the BAM
        for el in self.bam.pileup(contig_name, start, stop):
            pos=el.pos-start
            #### CHECK that POSITION is really inside the rage
            if pos>=0 and pos <= maxim :
                num_read=el.n
                data[0,pos]=el.pos
                data[1,pos]=num_read
                ### CHECK if the coverage is inside the BOUNDARIES
                if num_read>=min_coverage and num_read<max_coverage:
                    average_in_range+=num_read
                    sd_in_range+=num_read*num_read
                    hist_cov[int((num_read-min_coverage)/bin_size)]+=1
            else :
                bases_ouside_range+=1
        print("Bases outside RANGE:", bases_ouside_range)
        bases_in_range=np.sum(hist_cov)
        if normalized :
            hist_cov=hist_cov/(bases_in_range+0.0)
        sd_in_range=np.sqrt((sd_in_range - average_in_range*average_in_range/bases_in_range)/bases_in_range)
        average_in_range/=bases_in_range
        if self.verbose:
            print(f"Region {contig_name},{start}:{stop} [{min_coverage}:{max_coverage}] bases:{bases_in_range} avg:{average_in_range} stddev:{sd_in_range}")
        return hist_cov, bases_in_range, average_in_range, sd_in_range, data

    def validateCNV(self, contig, start, stop, min_cov=0, max_cov=500, bin_size=1, coverage_treshold=10., overlap_treshold=0.3, plot_hisogram=False, outputplot=None):
        lengthSV=stop-start
        ### Structural Variant
        sv=list(self.histogram(contig, start, stop, min_coverage=min_cov, max_coverage=max_cov, bin_size=bin_size, normalized=True))

        ### Surrounding Region
        L_ext=np.maximum(int(lengthSV/2),1000)
        ### Before - only positive position are accepted
        r1=list(self.histogram(contig, np.maximum(0,start-L_ext), start, min_coverage=min_cov, max_coverage=max_cov, bin_size=bin_size, normalized=False))
        ### After
        r2=list(self.histogram(contig, stop, stop+L_ext, min_coverage=min_cov, max_coverage=max_cov, bin_size=bin_size, normalized=False))
        ### Average in the Surrounding
        r0_base_count=r1[1]+r2[1]
        r0_avg=(r1[1]*r1[2]+r2[1]*r2[2])/(r1[1]+r2[1])
        r0_stddev=0.5*(r1[3]+r2[3])
        ### Normalize - Surrounding Area Histogram
        r_ext=(r1[0]+r2[0])/(r0_base_count+0.0)

        overlap=np.minimum(sv[0],r_ext)
        overlap_area=np.sum(overlap)

        check_coverage_Before=r1[2]>coverage_treshold
        check_coverage_After=r2[2]>coverage_treshold
        check_average_distance=abs(sv[2]-r0_avg)>r0_stddev
        check_overlapping_dist=overlap_area<overlap_treshold

        if self.verbose:
            print(f"CoverageBef[{r1[2]} > {coverage_treshold}]: {check_coverage_Before}")
            print(f"CoverageAft[{r2[2]} > {coverage_treshold}]: {check_coverage_After}")
            print(f"AverageDist[{r0_stddev} < |{sv[2]} - {r0_avg}|]: {check_average_distance}")
            print(f"Overlapping[{overlap_area} < {overlap_treshold}]: {check_overlapping_dist}")

        if plot_hisogram and r0_base_count+sv[1]>0:
            self.plotHistogram(sv, r1, r2, r_ext, overlap, bin_size, min_cov, max_cov)
            if outputplot is None :
                mp.show()
            else :
                mp.savefig(outputplot)
                mp.close()
                mp.cla()
                mp.clf()

        return check_coverage_Before and check_coverage_After and ( check_average_distance or check_overlapping_dist )

    def validateCNV_doublePlot(self, contig, start, stop, min_cov=0, max_cov=500, bin_size=1, coverage_treshold=10., overlap_treshold=0.3, plot_hisogram=False, outputplot=None):
        lengthSV=stop-start
        ### Structural Variant
        sv=list(self.histogram_and_coverage(contig, start, stop, min_coverage=min_cov, max_coverage=max_cov, bin_size=bin_size, normalized=True))

        ### Surrounding Region
        L_ext=np.maximum(int(lengthSV/2),1000)
        ### Before - only positive position are accepted
        r1=list(self.histogram_and_coverage(contig, np.maximum(0,start-L_ext), start, min_coverage=min_cov, max_coverage=max_cov, bin_size=bin_size, normalized=False))
        ### After
        r2=list(self.histogram_and_coverage(contig, stop, stop+L_ext, min_coverage=min_cov, max_coverage=max_cov, bin_size=bin_size, normalized=False))
        ### Average in the Surrounding
        r0_base_count=r1[1]+r2[1]
        r0_avg=(r1[1]*r1[2]+r2[1]*r2[2])/(r1[1]+r2[1])
        r0_stddev=0.5*(r1[3]+r2[3])
        ### Normalize - Surrounding Area Histogram
        r_ext=(r1[0]+r2[0])/(r0_base_count+0.0)

        overlap=np.minimum(sv[0],r_ext)
        overlap_area=np.sum(overlap)

        check_coverage_Before=r1[2]>coverage_treshold
        check_coverage_After=r2[2]>coverage_treshold
        check_average_distance=abs(sv[2]-r0_avg)>r0_stddev
        check_overlapping_dist=overlap_area<overlap_treshold

        if self.verbose:
            print(f"CoverageBef[{r1[2]} > {coverage_treshold}]: {check_coverage_Before}")
            print(f"CoverageAft[{r2[2]} > {coverage_treshold}]: {check_coverage_After}")
            print(f"AverageDist[{r0_stddev} < |{sv[2]} - {r0_avg}|]: {check_average_distance}")
            print(f"Overlapping[{overlap_area} < {overlap_treshold}]: {check_overlapping_dist}")

        if plot_hisogram and r0_base_count+sv[1]>0:
            self.plotHistogramAndCoverage(sv, r1, r2, r_ext, overlap, bin_size, min_cov, max_cov, contig=contig)
            if outputplot is None :
                mp.show()
            else :
                mp.savefig(outputplot)
                mp.close()
                mp.cla()
                mp.clf()

        return check_coverage_Before and check_coverage_After and ( check_average_distance or check_overlapping_dist )

    def plotHistogram(self, sv, r1, r2, r_ext, overlap, bin_size, min_cov, max_cov):
        ### Normalize
        r1[0] = r1[0] / r1[1]
        r2[0] = r2[0] / r2[1]
        max_cov_tmp=0
        for idx,val in enumerate(sv[0]+r1[0]+r2[0]) :
            if val>0 :
                max_cov_tmp=idx
        max_cov_tmp=min([(max_cov_tmp+1)*bin_size*1.25, max_cov])

        x = np.arange(min_cov, max_cov, bin_size)
        mp.plot(x, sv[0], '-', color='#ff0000', label='SV', linewidth=1)
        mp.plot(x, r1[0], '--', color='#0099ff', label='before', linewidth=0.5)
        mp.plot(x, r2[0], '--', color='#00cc99', label='after', linewidth=0.5)
        mp.plot(x, r_ext, '-', color='#00ff00', label='surrounding', linewidth=1)
        mp.plot(x, overlap, '-', color='#0000ff', linewidth=1)
        mp.xlim([min_cov,max_cov_tmp])
        mp.xlabel("Coverage")
        mp.ylabel("Num. Bases")

    def plotHistogramAndCoverage(self, sv, r1, r2, r_ext, overlap, bin_size, min_cov, max_cov, contig=None):
        min_indx = 0 #int(min_cov/bin_size)
        max_indx = int(np.ceil((max_cov - min_cov) / bin_size))

        pos_min=r1[4][0,0]
        pos_max=r2[4][0,-1]

        max_cov_tmp=min([max([max(sv[4][1,:]),max(r1[4][1,:]),max(r2[4][1,:])])*1.25, max_cov])
        ### Normalize
        r1[0] = r1[0] / r1[1]
        r2[0] = r2[0] / r2[1]
        max_hist_tmp=min(1., max([max(sv[0][min_indx:max_indx]),max(r1[0][min_indx:max_indx]),max(r2[0][min_indx:max_indx])])*1.25 )

        mp.figure(figsize=(6, 3))
        G = gridspec.GridSpec(1, 3)
        #################################################################################################
        ### COVERAGE PLOT
        axes_1 = mp.subplot(G[0, :2])
        axes_1.plot(sv[4][0,:].astype(int),sv[4][1,:], '-' , color='#ff0000', label='SV', linewidth=1)
        axes_1.plot(r1[4][0,:].astype(int),r1[4][1,:], '--', color='#0099ff', label='before', linewidth=1.)
        axes_1.plot(r2[4][0,:].astype(int),r2[4][1,:], '--', color='#00cc99', label='after', linewidth=1.)
        axes_1.set_xlabel('position')
        axes_1.set_xlim([pos_min,pos_max])
        axes_1.set_ylabel('coverage')
        axes_1.set_ylim([min_cov, max_cov_tmp])
        axes_1.grid(linestyle='--', linewidth=0.5)
        if contig is not None :
            axes_1.set_title(contig, position=(0.5, 0.5))
        #mp.subplots_adjust(rspace = .001)

        #################################################################################################
        ### HISTOGRAM
        axes_2 = mp.subplot(G[0,2])
        x = np.arange(min_cov, max_cov, bin_size)
        axes_2.plot(r1[0][min_indx:max_indx], x, '--', color='#0099ff', label='before', linewidth=1.)
        axes_2.plot(r2[0][min_indx:max_indx], x, '--', color='#00cc99', label='after', linewidth=1.)
        axes_2.plot(overlap[min_indx:max_indx], x, '-', color='#0000ff', linewidth=1)
        axes_2.plot(r_ext[min_indx:max_indx], x, '-', color='#0000ff', label='surrounding', linewidth=1)
        axes_2.plot(sv[0][min_indx:max_indx], x, '-', color='#ff0000', label='SV', linewidth=1)
        #axes_2.set_xlabel('Num. bases')
        #mp.subplots_adjust(rspace = .001)
        axes_2.set_xlim([0,max_hist_tmp])
        axes_2.set_ylim([min_cov, max_cov_tmp])
        axes_2.grid(linestyle='--', linewidth=0.5)
        #axes_2.set_ylabel('y')
        #if contig is not None :
        #    axes_2.set_title(contig)
        axes_2.set_xlabel("Probability")
        #axes_2.ylabel("Num. Bases")

        mp.subplots_adjust(bottom=0.1, right=0.95, top=0.95)
        mp.subplots_adjust(wspace=0.0, hspace=0.0)
        mp.tight_layout()
        mp.show()


if __name__ == '__main__':
    #### How to generate an INDEX for a BAM file:
    # samtools index bamfile.bam
    #### How to extract a subset of a BAM
    #  samtools view -b -h input.bam "Chr10:18000-45500" -o output.bam
    """
    source="/home/flanduzzi@iit.local/TestData/ShortBAM_P01_3B_F_16633.bam"
    bam=BAM_Manager(source)
    bam.verbose=True
    contig="chr22"
    start=50600000 #23100000
    stop=50600100 #23200000
    length=stop-start
    target_coverage=30
    bin_size=find_binsize(length,target_coverage)
    print("BIN SIZE:", bin_size)
    isValid=bam.validateCNV(contig, start, stop, bin_size=1, min_cov=0, max_cov=500, plot_hisogram=True)
    isValid=bam.validateCNV_doublePlot(contig, start, stop, bin_size=1, min_cov=0, max_cov=500, plot_hisogram=True)
    print("Structural Variant Quality Chek: ", isValid)
    """
    source="/home/flanduzzi@iit.local/TestData/SG_PD-7141_P_chr14:92358033-92458166.bam"
    bam=BAM_Manager(source)
    bam.verbose=True
    contig="chr14"
    start=92408033
    stop=92408166
    length=stop-start
    target_cov=100
    bin_size=find_binsize(length,target_cov)
    print("BIN SIZE:", bin_size)
    isValid=bam.validateCNV(contig, start, stop, bin_size=bin_size, min_cov=0, max_cov=500, plot_hisogram=True)
    isValid=bam.validateCNV_doublePlot(contig, start, stop, bin_size=1, min_cov=0, max_cov=500, plot_hisogram=True)
    print("Structural Variant Quality Chek: ", isValid)
