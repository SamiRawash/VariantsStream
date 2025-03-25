#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"get_vcf_field.h"
#include<getopt.h>
#include<ctype.h>


int main(int argc, char *argv[]){
    char *line = NULL;
    size_t len = 0;
    FILE *vcf = NULL;
    ssize_t read = 0;
    char *field;
    char c;
    int hflag = 0;
    int aflag = 0;
    int bflag = 0;
    char *sflag = NULL;
    char *copied_line;

    /* Help string */
    char help[] = "Usage: get_pass_variants [OPTION]... VCF_file\n  "
        "-h\tshow help options";


    while ((c = getopt (argc, argv, "h")) != -1)
        switch (c) {
            case 'h':
                hflag = 1;
                puts(help);
                return 1;
            case '?':
                if (isprint(optopt))
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
                return 1;
            default:
                abort();
        }

    argc -= optind;
    argv += optind;


    /* Check if there is a VCF file argument */
    if (argc < 1) {
        fprintf(stderr, "VCF file name not specified\n");
        return -1;
    }

    /*If the VCF file is not present */
    if (!(vcf = fopen(argv[0], "r"))) {
        fprintf(stderr, "Coud not find file: %s\n", argv[0]);
        return 1;
    }


    while ((read = getline(&line, &len, vcf) ) != -1) {
        if (line[0] != '#') {
            copied_line = (char *) malloc(strlen(line));
            strncpy(copied_line, line, strlen(line));
            field = get_vcf_field(copied_line, FILTER);
            if (field == NULL) {
                fprintf(stderr, "%s %s", argv[0], "is not a VCF file\n");
                exit(EXIT_FAILURE);
            }
            if (strcmp(field, "PASS") == 0) {
                fprintf(stdout, "%s", line);
            }
            free(copied_line);
        } else {
            fprintf(stdout, "%s", line);
        }

        
    }

    free(line);

    fclose(vcf);

    return 0;
}
