#pragma once

void graceful_shutdown(int exit_code);
void log_printf(char *formatString, ...);
int touch(char *filename);