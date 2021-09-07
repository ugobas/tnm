FILE *OpenFile(char *filename, char *envvar, char *mode, int *noenv);
char *GetWord(char *buffer, char *word, int maxlen, int comma);
void FreeArray2D(char   **array, int    dim1,  int    dim2);
char **Array2D(int size, int dim1,  int dim2);

#define KILLLEADSPACES(y,x)                                           \
                 do \
                 {  for((y)=(x); *(y) == ' ' || *(y) == '\t'; (y)++) ; } \
                 while(0)

#define TERMINATE(x) do {  int _terminate_macro_j;                    \
                        for(_terminate_macro_j=0;                     \
                            (x)[_terminate_macro_j];                  \
                            _terminate_macro_j++)                     \
                        {  if((x)[_terminate_macro_j] == '\n')        \
                           {  (x)[_terminate_macro_j] = '\0';         \
                              break;                                  \
                     }  }  }  while(0)
#define TOGGLE(x) (x) = (x) ? 0 : 1
