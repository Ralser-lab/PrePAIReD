#' Function used to import DIA-NN into R
#' @param diann_report_file Path to the DIA-NN report file
#' @param norm_intensity Boolean; whether to use raw or normalized intensities from the DIA-NN report
#' @param fasta_file Path to the FASTA files used to annotate the protein IDs to genes
#' @param organism Name of organism; used to extract protein description from the FASTA file
#' @param metadata_file Path to the metadata file
#' @return List object with peptide 

function import_dataset(diann_report_file, norm_intensity, fasta_file, organism = "HUMAN", metadata_file){
    require(msdap)
    # Import DIA-NN report file
    dataset <- import_dataset_diann(filename = diann_report_file, use_normalized_intensities = norm_intensity, do_plot = F)

    # Annotate proteins in the data with gene names
    dataset <- import_fasta(dataset, files = fasta_file)
    # Extract gene descriptions from the FASTA entries
    dataset$proteins$description <- sapply(dataset$proteins$fasta_headers, function(header){
        if(!grepl(";", header)){
            substring(header, first = gregexpr(pattern =paste("_", toupper(organism)),header)[[1]][1]+7, last = gregexpr(pattern ="OS",header)[[1]][1]-2)
        }else{
            paste(Reduce("c", lapply(strsplit(header, ";")[[1]], function(subheader) substring(subheader, first = gregexpr(pattern =paste("_", toupper(organism)),subheader)[[1]][1]+7, last = gregexpr(pattern ="OS",subheader)[[1]][1]-2))), collapse=";")
        }
    })
    
    # Import sample metadata
    dataset <- import_sample_metadata(dataset, filename = metadata_file)
    return(dataset)
}

function main(){
    #!/usr/bin/env Rscript
    library(optparse)
    
    option_list = list(
        make_option(c("-fd", "--file_diann"), type="character", default=NULL, help="Path to the DIA-NN report file", metavar="character"),
        make_option(c("-ni", "--norm_intensity"), type="logical", default=TRUE, help="logical; whether to use raw or normalized intensities from the DIA-NN report [default= %default]", metavar="logical")
    )
    
    opt_parser = OptionParser(option_list=option_list)
    opt = parse_args(opt_parser)
}
