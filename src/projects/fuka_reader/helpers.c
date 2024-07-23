#define _GNU_SOURCE
#include "bam.h"
#include "helpers.h"

void graceful_shutdown(int exit_code)
{
    fflush(stdout);
    bampi_finalize(0, 0);
    exit(exit_code);
}

void log_printf(char *formatString, ...)
{
    va_list args;
    va_start(args, formatString);
    char *message;
    vasprintf(&message, formatString, args);
    va_end(args);
    printf("fuka_reader: %s\n", message);
    free(message);
}

int touch(char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        return -1;
    }
    fclose(file);
    return 0;
}