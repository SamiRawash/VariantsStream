process CopyFollowingLinks {
	tag "MoveDuplsMetr"
input:
        path ff
        path destination
output:
        path "${destination}/${ff}"
script:
        """
        rsync -aLK ${ff} ${destination}
        """
}

process MultiQC_DuplMetr {
        publishDir "${multiqc_folder}", mode: 'move', overwrite: true 
input :
        path reports 
        path multiqc_folder
output :
        path "multiqc_report.html"
script:
        """
        multiqc .
        """
}