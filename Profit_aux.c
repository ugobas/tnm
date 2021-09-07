#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Profit_aux.h"

/*************************************************************************

   Program:
   File:       array.c, OpenFile.c, GetWord.c

   Version:    V1.4R
   Date:       18.03.94
   Function:   Allocate and free 2D arrays

   Copyright:  (c) SciTech Software 1993-4
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk

**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a
   commercial product except as described in the file COPYING.DOC.

**************************************************************************/

/************************************************************************/
/*>char *GetWord(char *buffer, char *word, int maxlen,  int    comma       Treat commas like white space?)
   -----------------------------------------------------------------
   Input:   char    *buffer     Input buffer to read words from
            int     maxlen      Max length of output word
	    int    comma       Treat commas like white space?
   Output:  char    *word       Word read from buffer
   Returns: char    *           Pointer to start of next word in buffer
                                or NULL

   Reads a whitespace delimted word out of buffer into word. If comma is
   1, then commas are treated just like white space, otherwise they
   are treated like normal characters.

   Words containing white space may be wrapped in double inverted commas.
   A \ is used as an escape character and maybe used to escape *any*
   following character. In particular:
      "\\" -> '\'     To get a backslash
      "\ " -> ' '     To get a hard whitespace (alternatively wrap the
                      string in double inverted commas)
      "\"" -> '"'     To get a double inverted comma

   10.06.99 Original   By: ACRM (based on code from Bioplib)
*/
char *GetWord(char *buffer, char *word, int maxlen, int comma)
{
   int  i, j;
   int dic    = 0,
        escape = 0;
   char *chp;

   /* Decrement maxlen so we can terminate correctly                    */
   maxlen--;

   /* Check validity of passed pointers                                 */
   if(word==NULL)
      return(NULL);

   word[0] = '\0';
   if(buffer==NULL)
      return(NULL);

   KILLLEADSPACES(chp, buffer);

   /* Run through each character in the input buffer                    */
   for(i=0, j=0; chp[i]; i++)
   {
      switch(chp[i])
      {
      case '\\':
         /* Use backslash as an escape character. If we've just had an
            escape, then simply store it
         */
         if(escape)
         {
            escape = 0;
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            escape = 1;
         }
         break;
      case '\"':
         /* Double inverted commas enclose strings containing white space
            If we've just had an escape then handle as a normal character,
            otherwise, toggle the dic flag
         */
         if(escape)
         {
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            TOGGLE(dic);
         }
         escape = 0;
         break;
      case ',':
         /* A comma is handled as white space or a normal character,
            depending on the comma flag
         */
         if(!comma)   /* Treat as default                               */
         {
            if(j<maxlen)
               word[j++] = chp[i];
            escape = 0;
            break;
         }
         /* Otherwise, if comma is true, just fall through to treat it
            like whitespace
         */
      case ' ':
      case '\t':
         /* If we are in double inverted commas or last char was an escape
            just handle as a normal character
         */
         if(dic || escape)
         {
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            /* Otherwise, this terminates the word, so terminate, move
               the pointer on and return
            */
            word[j] = '\0';
            chp += i;
            KILLLEADSPACES(chp, chp);
            if(comma)
            {
               /* If we are handling commas as whitespace, then k
                  the comma if found
               */
               if(*chp == ',') chp++;
            }
            if(*chp == '\0') chp = NULL;
            return(chp);
         }
         escape = 0;
         break;
      default:
         /* A normal character, copy it across                          */
         if(j<maxlen)
            word[j++] = chp[i];
         escape = 0;
      }
   }

   word[j] = '\0';
   return(NULL);
}

/************************************************************************/
/*>FILE *OpenFile(char *filename, char *envvar, char *mode, int *noenv)
   ---------------------------------------------------------------------
   Input:     char    *filename     Filename to be opened
              char    *envvar       Unix/MS-DOS environment variable
                                    Other OS assign name (with :)
              char    *mode         Mode in which to open file (r, w, etc)
   Output:    int    *noenv        Set to 1 under Unix/MS-DOS if
                                    the reason for failure was that the
                                    environment variable was not set.
   Returns:   FILE    *             File pointer or NULL on failure

   Attempts to open a filename as specified. Returns a file
   pointer. If this fails:

   Under UNIX/MS-DOS:
   gets the contents of the envvar environment variable and prepends
   that to the filename and tries again. If envvar was not set, noenv
   is set to 1 and the routine returns a NULL pointer.

   Under other OSs:
   prepends the envvar string onto the filename and tries to open the
   file again.

   Returns the pointer returned by the open() command after all this.

   22.09.94 Original    By: ACRM
   11.09.94 Puts a : in for the assign type.
   24.11.94 Added __unix define. Checks for trailing / in environment
            variable
   08.03.95 Corrected basename to filename in non-unix version
   09.03.95 Checks that filename is not a NULL or blank string
   28.07.05 Added conditionals for Mac OS/X: __MACH__ and __APPLE__
*/
FILE *OpenFile(char *filename, char *envvar, char *mode, int *noenv)
{
   char *datadir,
        buffer[160];
   FILE *fp;

   if(filename == NULL || filename[0] == '\0')
      return(NULL);

   if(noenv != NULL) *noenv = 0;

   /* Try to open the filename as specified                             */
   if((fp=fopen(filename,mode)) == NULL)
   {
      /* Failed, so build alternative directory/filename                */
#if (unix || __unix__ || MS_WINDOWS || __unix || __MACH__ || __APPLE__)
      if((datadir = getenv(envvar)) != NULL)
      {
         if(datadir[strlen(datadir)-1] == '/')
            sprintf(buffer,"%s%s",datadir,filename);
         else
            sprintf(buffer,"%s/%s",datadir,filename);
         fp = fopen(buffer,mode);
      }
      else
      {
         if(noenv != NULL) *noenv = 1;
         return(NULL);
      }
#else
      sprintf(buffer,"%s:%s",envvar,filename);
      fp = fopen(buffer,mode);
#endif
   }

   return(fp);
}

/************************************************************************/
/*>char **Array2D(int size, int dim1, int dim2)
   --------------------------------------------
   Input:   int   size    Size of an array element
            int   dim1    First dimension (number of rows)
            int   dim2    Second dimension (number of columns)
   Returns: char  **      Array of pointers. Must be cast to required
                          type

   Create a 2D array of elements of size `size' with dimensions `dim1'
   rows by `dim2' columns.

   07.10.92 Original
   12.07.93 Tidied and commented
*/
char **Array2D(int size,
               int dim1,
               int dim2)
{
   char  **array  = NULL;
   int   i;

   /* Allocate memory for the outer dimension array                     */
   if((array = (char **)malloc(dim1 * sizeof(char *))) == NULL)
      return(NULL);

   /* Set all positions to NULL                                         */
   for(i=0; i<dim1; i++)   array[i] = NULL;

   /* Allocate memory for each array in the second dimension            */
   for(i=0; i<dim1; i++)
   {
      /* If allocation fails, jump to badexit                           */
      if((array[i] = (char *)malloc(dim2 * size)) == NULL)
         goto badexit;
   }

   return(array);

badexit:
   for(i=0; i<dim1; i++)   if(array[i]) free(array[i]);
   free(array);
   return(NULL);
}

/************************************************************************/
/*>void FreeArray2D(char **array, int dim1, int dim2)
   --------------------------------------------------
   Input:   char  **    Array of pointers to be freed
            int   dim1  First dimension (number of rows)
            int   dim2  Second dimension (number of columns)

   Frees a 2D array with dimensions `dim1' rows by `dim2' columns.

   07.10.92 Original
*/
void FreeArray2D(char   **array,
                 int    dim1,
                 int    dim2)
{
   int   i;

   if(array)
   {
      for(i=0; i<dim1; i++)   if(array[i]) free(array[i]);
      free(array);
   }
}
