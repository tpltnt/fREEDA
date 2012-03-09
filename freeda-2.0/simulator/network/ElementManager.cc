#include "Circuit.h"
#include "../network/ADInterface.h"
#include "../elements/element_headers.h"
#include "../inout/environment.h"

extern const char* freeda_version;
extern void sanityCheck();

ElementManager* ElementManager::em = NULL;

ElementManager::ElementManager()
{
  // Do nothing
}

ElementManager::~ElementManager()
{
  // Do nothing
}

// Only one element manager can exits.
ElementManager* ElementManager::getElementManager()
{
  if (!em)
    em = new ElementManager;
  return em;
}

Element* ElementManager::createElement(const string& elem_type,
const string& iname)
throw(string&)
{
  //  string instance_name = elem_type + ":" + iname;
  Element* new_elem = NULL;
  // Find the type of element (code automatically generated).
	#include "../elements/create_element.cc"

  if (new_elem)
    return new_elem;
  else
    throw("Element type " + elem_type + " does not exist.");
}


void ElementManager::printCatalog() const
{
  char pdfFileRoot[FILENAME_MAX];
  char htmlFilename[FILENAME_MAX];
  char elementHtmlFilename[FILENAME_MAX];
  char elementName[FILENAME_MAX];
  int useLocalElementHtmlDocumentation;
  int useLocalElementPdfDocumentation;
  int length;
  FILE *test_F;


  /* Output environment variables */
  printf("**  Environment variables: **\n");
  printf("FREEDA_HOME = %s\n",env_freeda_home);
  printf("FREEDA_LIBRARY = %s\n",env_freeda_library);
  printf("FREEDA_PROJECTS = %s\n",env_freeda_projects);
  printf("FREEDA_PATH = %s\n",env_freeda_path);
  printf("FREEDA_BIN = %s\n",env_freeda_bin);
  printf("FREEDA_SIMULATOR = %s\n",env_freeda_simulator);
  printf("FREEDA_ELEMENTS = %s\n",env_freeda_elements);
  printf("FREEDA_DOCUMENTATION = %s\n",env_freeda_documentation);
  printf("FREEDA_WEB_DOCUMENTATION = %s\n",env_freeda_web_documentation);
  printf("FREEDA_BROWSER = %s\n",env_freeda_browser);
  printf("\n");

  // Vector to hold the original object for each element type.
  Element** elem_vector = NULL;
  // Allocate memory for the vector
  elem_vector = new Element*[ELEM_TYPES];
  // Allocate one object each element type
  // The code for this is automatically generated
  #include "../elements/create_dummy_elem.cc"

  // Set up documentation for individual elements
  {
  char s[6*FILENAME_MAX];
  // Check to see if env_freeda_web_documentation is local.
  sprintf(s,"test -d %s",env_freeda_web_documentation);
  if(!system(s)) // It is local.
    useLocalElementPdfDocumentation=1;
  else // It is probably a url.
    useLocalElementPdfDocumentation=0;

  // Check to see if directory for element documentation exists.
  sprintf(s,"test -d %s/elements",env_freeda_documentation);
  if(system(s))
    { // Directory does not exits, so create it.
    sprintf(s,"mkdir -p %s/elements",env_freeda_documentation);
    if(system(s))
      { // Cannot create directory
      printf("Cannot create %s/elements: cannot create local documentation \n",
        env_freeda_documentation);
      printf("This is a fatal error.");
      exit(1);
      }
    }
  useLocalElementHtmlDocumentation=1;


  if(useLocalElementPdfDocumentation)
    {
    // copy pdfs from elements tree to documentation.
    // Check to see if source tree exists.
    sprintf(s,"test -d %s",env_freeda_elements);
    if(system(s))
      {
      printf("Local pdf's not found. Cannot open: %s\n",env_freeda_elements);
      printf("Using web for pdf documentation instead.\n");
      useLocalElementPdfDocumentation=0;
      }
    else
      {
      printf("Copying pdf files from the tree: %s\n",env_freeda_elements);
      printf("  and putting them in: %s\n",env_freeda_web_documentation);
      sprintf(s,"find %s -name \"*.pdf\" -exec cp {} %s/elements \\;",
        env_freeda_elements, env_freeda_web_documentation);
      if(system(s))
        {
        printf("Fatal error: operation to gather pdf files failed.\n");
        exit(1);
        }
      }
    }

  // Now point to the location of the element pdfs.
  sprintf(pdfFileRoot,"%s/elements",env_freeda_web_documentation);
  printf("Root for PDF documentation for individual elements ");
  printf("is\n  %s .\n",pdfFileRoot);
  }

  // Open catalog html file
  // Use documentation environment variable
  sprintf(htmlFilename,"%s/fr_elements.html",env_freeda_documentation);
  FILE *out_f = fopen(htmlFilename, "w");
  if (out_f)
    {
    printf("** Creating catalog in %s\n",htmlFilename);
    }
  else
    {
    fprintf(stderr, "Could not open %s for writing.\n", htmlFilename);
    exit(1);
    }

  fprintf(out_f,
    "<HTML> \n <HEAD><TITLE>fREEDA(TM) Element Catalog</TITLE></HEAD> \n");
  fprintf(out_f,
    "<HR SIZE=4>\n<h3 ALIGN=LEFT>fREEDA(TM) Element Catalog</h3>\n\n");
  fprintf(out_f,
    "<p ALIGN=LEFT>Version %s compiled on %s   %s</p>\n<HR SIZE=4>\n",
    freeda_version,__DATE__,__TIME__);
  fprintf(out_f, "<table width=\"100%%\" border=\"0\">\n");
  fprintf(out_f, "<tr>\n");
  fprintf(out_f,
    "<th width=\"20%%\" scope=\"col\"><p align=\"left\">ELEMENT</p></th>\n");
  fprintf(out_f,
   "<th width=\"80%%\" scope=\"col\"><p align=\"left\">DESCRIPTION</p></th>\n");
  fprintf(out_f, "</tr>\n");

  for (int i = 0; i < ELEM_TYPES; i++)
    {
    Element* elem = elem_vector[i];
    elementName[0]=0;
    strcpy(elementName,elem->getName().c_str());
    fprintf(out_f, "<tr>\n<td>"); //Start the row entry
    // hyperref element's html filename.
    sprintf(elementHtmlFilename,"elements/%s.html", elementName);
    fprintf(out_f,"<a href=\"%s\">",elementHtmlFilename);
    fprintf(out_f, "%s</a></td>\n", elementName);
    fprintf(out_f, "<td>%s</td>\n</tr>\n", elem->getDescription().c_str());
    printCatalogElement(elementName,0); // Create element's html file.
    }

  fprintf(out_f, "</table>\n");

  fprintf(out_f, "</BODY>\n</HTML>\n");
  fclose(out_f);

  {
  char s[FILENAME_MAX];
  printf("Opening %s\n",htmlFilename);
  sprintf(s,"%s %s",env_freeda_browser,htmlFilename);
  system(s);
  }

  delete [] elem_vector; // Free the dummy elements
}

void ElementManager::printCatalogElement(char* elementName,
  int launchBrowserFlag) const
{
  char htmlFilename[FILENAME_MAX];
  char pdfFilename[FILENAME_MAX];
  char fullpdfFilename[FILENAME_MAX];
  char thisElementName[FILENAME_MAX];
  char s[3*FILENAME_MAX];
  Element* elem;
  int elementIndex;
  int length;
  int localPdfFile;
  FILE *test_F;
  FILE *out_f;
  Element** elem_vector = NULL;
  elem_vector = new Element*[ELEM_TYPES];
  #include "../elements/create_dummy_elem.cc"

  // Get element
  for (int i = 0; i < ELEM_TYPES; i++)
    {
    elem = elem_vector[i];
    thisElementName[0]=0;
    strcpy(thisElementName,elem->getName().c_str());
    if(!strcmp(elementName,thisElementName))
      {
      elementIndex = i;

      // open output file
      {
      // Use documentation environment variable
      length = strlen(env_freeda_documentation) + 21 + strlen(thisElementName);
      if(length >= FILENAME_MAX)
        {
        fprintf(stderr," ** Filename too long.\n");
        return;
        }

      // Check to see if directory for element documentation exists.
      sprintf(s,"test -d %s/elements",env_freeda_documentation);
      if(system(s))
        { // Directory does not exits, so create it.
        sprintf(s,"mkdir -p %s/elements",env_freeda_documentation);
        if(system(s))
          { // Cannot create directory
          printf("Cannot create %s/elements: \n",env_freeda_documentation);
          printf("This is a fatal error.");
          exit(1);
          }
        }


      htmlFilename[0] = 0;
      strcpy(htmlFilename,env_freeda_documentation);
      strcat(htmlFilename,"/elements/");
      strcat(htmlFilename,thisElementName);
      strcat(htmlFilename,".html");
      out_f = fopen(htmlFilename, "w");
      if (out_f == NULL)
        {
        fprintf(stderr, "Could not open %s for writing.\n", htmlFilename);
        return;
        }
      else
        {
        printf("** Creating %s\n",htmlFilename);
        }
      }

       fprintf(out_f,
          "<HTML> \n <HEAD><TITLE>%s fREEDA(TM) Element</TITLE></HEAD> \n",
          thisElementName);

       fprintf(out_f,"<HR SIZE=4>\n");

       fprintf(out_f, "<table width=\"100%%\" border=\"0\">\n");
       fprintf(out_f, "<tr>\n");
       fprintf(out_f,
         "<th width=\"20%%\" scope=\"col\"><h2 ALIGN=LEFT>%s</h2></th>\n",
         thisElementName);
       fprintf(out_f,
         "<th width=\"80%%\" scope=\"col\"><p ALIGN=RIGHT>%s</p></th>\n",
         elem->getDescription().c_str());
       fprintf(out_f, "</tr>\n");
       fprintf(out_f, "</table>\n");
       fprintf(out_f,"<HR SIZE=4>\n");
       fprintf(out_f, "<p ALIGN=LEFT>\n");


       // see if pdf documentation exists and hyperref
       {
       sprintf(pdfFilename,"%s.pdf",thisElementName);
       // check to see if it exists locally
       sprintf(fullpdfFilename,"%s/elements/%s",
         env_freeda_documentation ,pdfFilename);
       test_F = fopen(fullpdfFilename, "r");
       if(test_F)
         {
         // The documentation is available locally.
         fprintf(out_f,"<a href=\"%s\">",pdfFilename);
         }
       else
         {
         // point to the on-line documentation
         sprintf(fullpdfFilename,"%s/elements/%s",
           env_freeda_web_documentation, pdfFilename);
         fprintf(out_f,"<a href=\"%s\">",fullpdfFilename);
         }
       fclose(test_F);
       }


       fprintf(out_f, "Click here for pdf documentation.</a></p>\n");

       fprintf(out_f, "Author(s): %s<br>\n", elem->getAuthor().c_str());

       if (elem->satisfies(MULTI_REF))
         fprintf(out_f, "<li>Multi-referenced element. </li><br>\n");
       if (!elem->getNumTerms())
         fprintf(out_f,
           "<li>Number of terminals for this element is variable. </li><br>\n");

       fprintf(out_f, "<br>\n<B>Usage:</B> <br>\n");
       fprintf(out_f,
         "<ul><B>%s</B>:&lt;instance name&gt; ",elem->getName().c_str());

       if (elem->getNumTerms())
         {
         for (unsigned j = 0; j < elem->getNumTerms(); j++)
           {
           fprintf(out_f, "<B>n%d</B> ", j+1);
           }
         }
       else
         {
         fprintf(out_f, "<B>n1 n2 ... </B>");
         }

       // Special case for subcircuit instance
       if (elem->getName() == "x")
         fprintf(out_f, "&lt;subcircuit name&gt; ");

       if (elem->getNumberOfParams())
         {
         fprintf(out_f, " &lt;parameter list&gt; </ul><br>\n");
         fprintf(out_f, "<table border cellpadding=5 cellspacing=0> \n\n");
         fprintf(out_f, "<tr> <th>Parameter</th> <th>Type</th>");
         fprintf(out_f, "<th>Default value</th> <th>Required?</th> </tr> \n");

         for (unsigned j = 0; j < elem->getNumberOfParams(); j++)
           {
           string parname, comment, type, dflt_val, required;
           elem->getParamDesc(j, parname, comment, type, dflt_val, required);
           fprintf(out_f, "<tr><td><B>%s:</B> %s</td>\n",
           parname.c_str(), comment.c_str());;
           fprintf(out_f, "<td> %s</td>\n", type.c_str());;
           fprintf(out_f, "<td> %s</td>\n", dflt_val.c_str());
           fprintf(out_f, "<td> %s</td></tr>\n\n", required.c_str());
           }
         fprintf(out_f, "</table>\n");
         }
       else
         fprintf(out_f, "</ul><br><br>\n");

      fprintf(out_f,
        "<p ALIGN=LEFT>fREEDA Version %s compiled on %s  %s</p>\n<HR SIZE=4>\n",
        freeda_version,__DATE__,__TIME__);

      fclose(out_f);
      break;
      }
    }

  if(launchBrowserFlag)
    {
    length = strlen(env_freeda_browser) +1 + strlen(htmlFilename) +1;
    if(length >= FILENAME_MAX)
      {
      printf(" ** Filename too long.\n");
      return;
      }
    s[0] = 0;
    strcpy(s,env_freeda_browser);
    strcat(s," ");
    strcat(s,htmlFilename);
    system(s);
    }

  // Free the dummy elements.  Only the row dimension was allocated.
  delete [] elem_vector;
}

