GeDi_footer <- fluidRow(
  column(
    width = 1,
    align = "right",
    a(
      href = "https://github.com/AnnekathrinSilvia/GeDi",
      target = "_blank",
      img(src = "GeDi/GeneTonic.png", height = "50px")
    )
  ),
  column(
    width = 11,
    align = "center",
    "GeDi is a project developed by Annekathrin Ludt
     in the Bioinformatics division of the ",
    tags$a(href = "http://www.unimedizin-mainz.de/imbei", "IMBEI"),
    "- Institute for Medical Biostatistics, Epidemiology and Informatics", br(),
    "License: ", tags$a(href = "https://opensource.org/licenses/MIT", "MIT"),
    "- The GeDi package is developed and available on ",
    tags$a(href = "https://github.com/AnnekathrinSilvia/GeDi", "GitHub")
  )
)
