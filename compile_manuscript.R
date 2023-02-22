compile_manuscript <- function() 
{
  tinytex::pdflatex("manuscript.tex")
}

open_pdf <- function(file) 
{
  if (exists("RStudio.Version") && RStudio.Version()$mode != "server")
  {
    switch(
      Sys.info()[['sysname']],
      Windows = {shell.exec(normalizePath(file))},
      Linux   = {system2(getOption("pdfviewer", default=''), args = file)},
      Darwin = {system2('open', args = file, wait = FALSE)}
    )
  }
}
