/* main.h */
/* Bernd Bruegmann, 12/99 */

void read_command_line(int argc, char **argv);
void parse_parameter_file(char *parfile);
int iterate_parameters(void);
void initialize_libraries(void);
void initialize_grid(tG *g);
void evolve_grid(tG *g);
void advance_levelandsublevels(tG *g, int l);
void advance_level(tG *g, int l);
void handle_ghost_level(tG *g, int iteration);
void analyze_level(tG *g, int l);
void solve_elliptic(tG *g, int l);
void finalize_grid(tG *g);
void finalize_libraries(void);
void free_stuff(tG *g);


/* physics.c */
void init_physical_objects(tG* grid);






void freeFuns();
void freeParameters();
void freeparameters();
void freeVars();


