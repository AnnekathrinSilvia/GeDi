#' GeDi main function
#'
#' @param genesets a `data.frame`, The input data used for GeDi. This should be
#'                 a `data.frame` of at least two columns.One column should be
#'                 called "Genesets" and contain some sort of identifiers for
#'                 the individual genesets. In this application, we use the term
#'                 "Genesets" to refer to collections of individual genes, which
#'                 share common biological characteristics or functions. Such
#'                 genesets can for example be obtained from databases such as
#'                 the Gene Ontology (GO), the Kyoto Encyclopedia of Genes and
#'                 Genomes (KEGG), Reactome, or the Molecular Signatures
#'                 Database (MSigDB). The identifiers used in these databases
#'                 can be directly used as geneset identifiers in GeDi. The
#'                 second column should be called "Genes" and contain a list of
#'                 genes belonging to the individual genesets in the "Genesets"
#'                 column. In order to leverage all of the functionality
#'                 available in GeDi, the column has to contain gene names and
#'                 no other commonly used identifiers. The column names are case
#'                 sensitive.
#' @param ppi_df a `data.frame`, Protein-protein interaction (PPI) network data
#'               frame. The object is expected to have three columns, `Gene1`
#'               and `Gene2` which specify the gene names of the interacting
#'               proteins in no particular order (symmetric interaction) and a
#'               column `combined_score` which is a numerical value of the
#'               strength of the interaction.
#' @param distance_scores A [Matrix::Matrix()] of (distance) scores
#' @param gtl A `GeneTonicList`object generated with 
#'            [GeneTonic::GeneTonic_list()], containing the functional enrichment
#'            results.
#' @param col_name_genesets character, the name of the column in which the
#'                          geneset ids are listed. Defaults to "Genesets".
#' @param col_name_genes character, the name of the column in which the genes
#'                       are listed. Defaults to "Genes".
#'
#' @return A Shiny app object is returned
#' @export
#' @import visNetwork
#' @import shiny
#' @import shinyBS
#' @import fontawesome
#' @importFrom bs4Dash box bs4DashPage bs4DashNavbar bs4DashBrand bs4DashSidebar
#' bs4SidebarMenu bs4SidebarMenuItem bs4DashBody bs4DashControlbar bs4DashFooter
#' bs4TabItems bs4TabItem renderbs4InfoBox updateBox bs4Card
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom rintrojs introjs
#' @importFrom utils read.delim data
#' @importFrom shinycssloaders withSpinner
#' @importFrom igraph V degree delete_vertices get.edgelist
#'
#' @examples
#' if (interactive()) {
#'   GeDi()
#' }
#' # Alternatively, you can also start the application with your data directly
#' # loaded.
#'
#' data("macrophage_topGO_example", package = "GeDi")
#' if (interactive()) {
#'   GeDi(genesets = macrophage_topGO_example)
#' }
GeDi <- function(genesets = NULL,
                 ppi_df = NULL,
                 distance_scores = NULL,
                 gtl = NULL, 
                 col_name_genesets = "Genesets",
                 col_name_genes = "Genes") {
  oopt <- options(spinner.type = 6, spinner.color = "#0092AC")
  on.exit(options(oopt))

  usage_mode <- "shiny_mode"

  if (!(is.null(genesets))) {
    genesets <- .checkGenesets(genesets,
                               col_name_genesets,
                               col_name_genes)
  }
  if (!(is.null(ppi_df))) {
    ppi <- .checkPPI(ppi_df)
  }
  
  if (!(is.null(distance_scores))) {
    stopifnot("When providing distance scores, you also need to provide the geneset data" = !is.null(genesets))
    distance_scores <- .checkScores(genesets, distance_scores)
  }
  
  if(!(is.null(gtl))){
    genesets <- .checkGTL(gtl)
  }

  # UI definition -----------------------------------------------------------

  # dashpage definition -----------------------------------------------------
  gedi_ui <- bs4DashPage(
    title = "GeDi",
    dark = NULL,
    help = NULL,
    # navbar definition -------------------------------------------------------
    header = bs4DashNavbar(
      tagList(tags$code(tags$h3("GeDi"))),
      actionButton(
        "bookmarker",
        label = "Bookmark",
        icon = icon("heart"),
        style = .actionButtonStyle,
        class = "ml-5"
      ),
      tags$span(style = "display:inline-block; width: 2%"),
      tagList(
        shinyWidgets::dropdownButton(
          inputId = "ddbtn_info",
          circle = FALSE,
          status = "info",
          icon = icon("info"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(
            title = "More info on GeDi and on the current session"),
          tags$h5("Additional information"),
          actionButton(
            inputId = "btn_docs_vignette",
            icon = icon("book-open"),
            label = "Open GeDi Vignette",
            style = .actionButtonStyle,
            onclick = ifelse(
              system.file("doc", "GeDi_manual.html", package = "GeDi") != "",
              "",
              "window.open('https://annekathrinsilvia.github.io/GeDi/articles/GeDi_manual.html', '_blank')"
            )
          ),
          actionButton(
            inputId = "btn_info_session",
            icon = icon("circle-info"),
            label = "About this session",
            style = .actionButtonStyle
          ),
          actionButton(
            inputId = "btn_info_gedi",
            icon = icon("heart"),
            label = "About GeDi",
            style = .actionButtonStyle
          )
        )
      ),
      title = bs4DashBrand(title = HTML("<small>GeDi</small>"),
                           href = "https://github.com/AnnekathrinSilvia/GeDi",),
      skin = "light",
      border = FALSE,
      controlbarIcon = icon("gears"),
      fixed = TRUE
    ),

    # sidebar definition ------------------------------------------------------
    sidebar = bs4DashSidebar(
      id = "sidebar",
      title = HTML("<small>GeDi</small>"),
      src = "GeDi/GeDi.png",
      skin = "light",
      status = "info",
      brandColor = NULL,
      url = "https://github.com/AnnekathrinSilvia/GeDi",
      collapsed = TRUE,
      elevation = 1,
      opacity = 0.8,
      bs4SidebarMenu(
        id = "tabs",
        bs4SidebarMenuItem("Welcome!",
                           tabName = "tab_welcome",
                           icon = icon("house")),
        bs4SidebarMenuItem(
          "Data Input",
          tabName = "tab_data_input",
          icon = icon("database")
        ),
        bs4SidebarMenuItem(
          "Distance Scores",
          tabName = "tab_scores",
          icon = icon("table")
        ),
        bs4SidebarMenuItem(
          "Clustering Graph",
          tabName = "tab_graph",
          icon = icon("circle-nodes")
        ),
        bs4SidebarMenuItem(
          "Report",
          tabName = "tab_report",
          icon = icon("file-lines")
        )
      )
    ),

    # body definition ---------------------------------------------------------
    body = bs4DashBody(
      rintrojs::introjsUI(),
      tags$head(tags$style(
        HTML(
          ".shiny-output-error-validation {
            font-size: 15px;
            color: forestgreen;
            text-align: center;
            }
            #myScrollBox{
              overflow-y: scroll;
              .dataTables_wrapper{
                overflow-x: scroll;
              }
            }
          .tooltip > .tooltip-inner {
                 width: 1000px;
                 color: white;
                 background-color: black;
                 }
          .introjs-tooltip {
            max-width: 700px;
            min-width: 300px;
          }"
        )
      )),
      tags$head(
        tags$style(
          ".biocdlbutton{background-color:#0092AC;} .biocdlbutton{color: #ffffff;}"
        )
      ),
      tags$head(tags$style("#shiny-modal img { max-width: 100%; }")),
      tags$script(
        HTML(
          "$(function(){
      $(document).keyup(function(e) {
      if (e.which == 17) {
        $('#bookmarker').click()
      }
      });
      })"
        )
      ),
      bs4TabItems(
        # ui panel welcome ---------------------------------------------------
        bs4TabItem(tabName = "tab_welcome",
                   tagList(
                     fluidRow(column(width = 11),
                              column(
                                width = 1,
                                actionButton(
                                  "tour_welcome",
                                  label = "",
                                  icon = icon("circle-question"),
                                  style = .actionButtonStyle
                                )
                              )),
                     uiOutput("ui_panel_welcome")
                   )),
        # ui panel data input -----------------------------------------------
        bs4TabItem(
          tabName = "tab_data_input",
          tagList(
            fluidRow(
              column(width = 11),
              column(
                width = 1,
                br(),
                actionButton(
                  "tour_data_input",
                  label = "",
                  icon = icon("circle-question"),
                  style = .actionButtonStyle
                ),
                shinyBS::bsTooltip(
                  id = "tour_data_input",
                  title = "Click me to start a tour of this section!",
                  placement = "top",
                  trigger = "hover",
                  options = list(container = "body")
                ),
                style = "float:right"
              )
            ),
            uiOutput("ui_panel_data_input"),
            uiOutput("ui_filter_data"),
            uiOutput("ui_panel_specify_species"),
            uiOutput("ui_panel_download_ppi")
          )
        ),
        # ui panel scores -------------------------------------------------
        bs4TabItem(tabName = "tab_scores",
                   tagList(
                     fluidRow(
                       column(width = 11),
                       column(
                         width = 1,
                         br(),
                         actionButton(
                           "tour_scoring",
                           label = "",
                           icon = icon("circle-question"),
                           style = .actionButtonStyle
                         ),
                         shinyBS::bsTooltip(
                           id = "tour_scoring",
                           title = "Click me to start a tour of this section!",
                           placement = "top",
                           trigger = "hover",
                           options = list(container = "body")
                         ),
                         style = "float:right"
                       )
                     ),
                     uiOutput("ui_panel_scores")
                   )),
        # ui panel graph ---------------------------------------------------
        bs4TabItem(tabName = "tab_graph",
                   tagList(
                     fluidRow(
                       column(width = 11),
                       column(
                         width = 1,
                         br(),
                         actionButton(
                           "tour_graph",
                           label = "",
                           icon = icon("circle-question"),
                           style = .actionButtonStyle
                         ),
                         shinyBS::bsTooltip(
                           id = "tour_graph",
                           title = "Click me to start a tour of this section!",
                           placement = "top",
                           trigger = "hover",
                           options = list(container = "body")
                         ),
                         style = "float:right"
                       )
                     ),
                     uiOutput("ui_panel_graph")
                   )),
        bs4TabItem(tabName = "tab_report",
                   uiOutput("ui_panel_report"))
      )
    ),
    # footer definition -------------------------------------------------------
    footer = bs4DashFooter(left =
                             fluidRow(
                               column(
                                 width = 1,
                                 align = "right",
                                 a(
                                   href = "https://github.com/AnnekathrinSilvia/GeDi",
                                   target = "_blank",
                                   img(src = "GeDi/GeDi.png", height = "50px")
                                 )
                               ),
                               column(
                                 width = 11,
                                 align = "center",
                                 "GeDi is a project developed by Annekathrin
                                 Silvia Nedwed in the Bioinformatics division
                                 of the ",
                                 tags$a(
                                   href = "http://www.unimedizin-mainz.de/imbei",
                                   "IMBEI"),
                                 "- Institute for Medical Biostatistics,
                                 Epidemiology and Informatics",
                                 br(),
                                 "License: ",
                                 tags$a(
                                   href = "https://opensource.org/licenses/MIT",
                                   "MIT"),
                                 "- The GeDi package is developed and available
                                 on ",
                                 tags$a(
                                   href = "https://github.com/AnnekathrinSilvia/GeDi",
                                   "GitHub")
                               )
                             ),
                           right = NULL)
  )
  # nocov start
  # server definition ---------------------------------------------------------
  gedi_server <- function(input, output, session) {
    # initializing reactives --------------------------------------------------
    reactive_values <- reactiveValues()

    if (!is.null(genesets)) {
      reactive_values$genesets <- genesets
      reactive_values$gs_names <- genesets[[col_name_genesets]]
      reactive_values$genes <- getGenes(genesets, gene_name = col_name_genes)
      reactive_values$gs_description <-
        .getGenesetDescriptions(genesets)
    } else {
      reactive_values$genesets <- NULL
      reactive_values$gs_names <- NULL
      reactive_values$gs_description <- NULL
      reactive_values$genes <- NULL
    }

    if (!(is.null(ppi_df))) {
      reactive_values$ppi <- ppi_df
    } else {
      reactive_values$ppi <- NULL
    }

    reactive_values$filtered_genesets <- NULL
    reactive_values$alt_names <- FALSE
    reactive_values$species <- NULL

    if (!(is.null(distance_scores))) {
      reactive_values$scores[["My_Score"]] <- distance_scores
    } else {
      reactive_values$scores <- list()
    }
    reactive_values$seeds <- NULL
    reactive_values$cluster <- NULL

    # bookmarks
    reactive_values$bookmarked_genesets <- NULL
    reactive_values$bookmarked_cluster <- NULL

    # panel Welcome ----------------------------------------------------------
    output$ui_panel_welcome <- renderUI({
      tagList(fluidRow(column(
        width = 12,
        includeMarkdown(system.file("extdata", "welcome.md", package = "GeDi")),
        br(),
        uiOutput("ui_summary_data_generation")
      )),
      fluidRow(
        column(
          width = 12,
          includeMarkdown(
            system.file("extdata", "welcome_introduction.md", package = "GeDi")
          ),
          actionButton(
            inputId = "example_button",
            label = "Click me!",
            style = .actionButtonStyle
          ),
          br(),
          p(
            "you can click on it and it will start a process in the app.
          We tried to choose the button text wisely so that you can easily
          understand the process that will be started."
          ),
          p(
            "Besides buttons, you will also find input fields where you can type
          something like this"
          ),
          textInput(inputId = "exampleText",
                    label = "Write some text here"),
          p(
            "In these boxes you generally have to provide some information.
          In most cases, you will see that once you click inside the box, a
          drop-down menu will appear with suggestions for the specific field to
          provide you with an idea of what is expected."
          ),
          p(
            "Otherwise, there are also options, where you have to select something"
          ),
          radioButtons(
            inputId = "exampleRadioButton",
            label = "Select the option of your choice",
            choices = c("Option 1", "Option 2", "Option 3")
          ),
          p(
            "and sliders where you can define a certain value based on an
          interval of available values"
          ),
          sliderInput(
            inputId = "exampleSlider",
            label = "Choose a value",
            max = 100,
            min = 0,
            value = 50
          ),
          p(
            "We tried to keep the app as simple and self-explanatory as possible.
          However, if you happen to be stuck at some point, we provide for each
          panel of this app a dedicated tour which will explain each element
          of the respective panel in detail. You can access the tours via the
          question mark in the upper right corner, which you will find in each
          panel."
          )
        )
      ))
    })


    output$ui_summary_data_generation <- renderUI({
      box(
        width = 12,
        title = "Summary on data generation",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        htmltools::includeMarkdown(
          system.file("extdata", "summary_data_generation.md", package = "GeDi")
        )
      )
    })

    # panel Data Input ------------------------------------------------------
    output$ui_panel_data_input <- renderUI({
      tagList(
        box(
          width = 12,
          title = "Step 1",
          status = "danger",
          solidHeader = TRUE,
          h2("Provide your Genesets as input data"),
          fluidRow(
            column(
              width = 6,
              fileInput(
                inputId = "inputgenesetfile",
                label = "Input a geneset file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/tab-separated-values",
                  "text/plain",
                  ".RDS",
                  ".xslx"
                ),
                multiple = FALSE
              ),
              br(),
              "... or you can also ",
              br(),
              actionButton(
                "btn_loaddemo",
                "Load the demo data",
                icon = icon("play-circle"),
                class = "btn btn-info",
                style = .actionButtonStyle
              ),
              br(),
              "..and here you can have a look at how the data should be structured",
              actionButton(
                "btn_showDataStructure",
                "Have a look at the data structure",
                icon = icon("magnifying-glass"),
                class = "btn btn-info",
                style = .actionButtonStyle
              ),
              br(),
              p(),
              uiOutput("ui_alt_column_names")
            ),
            column(
              width = 6,
              box(
                id = "Genesets_preview",
                width = NULL,
                title = "Genesets preview",
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                collapsed = TRUE,
                fluidRow(column(
                  width = 12,
                  offset = 0.5,
                  DT::dataTableOutput("dt_genesets")
                ))
              )
            )
          )
        )
      )
    })

    readGenesetsTxt <- reactive({
      if (is.null(input$inputgenesetfile)) {
        return(NULL)
      }

      guessed_sep <-
        .sepguesser(input$inputgenesetfile$datapath)
      genesets <-
        utils::read.delim(
          input$inputgenesetfile$datapath,
          header = TRUE,
          as.is = TRUE,
          sep = guessed_sep,
          quote = "",
          check.names = FALSE
        )
      return(genesets)
    })

    output$ui_alt_column_names <- renderUI({
      if (reactive_values$alt_names == FALSE) {
        return(NULL)
      }
      fluidRow(
        column(
          width = 12,
          "It seems like your data does not contain columns named Genesets and
          Genes. Please select the column which contains the identifiers for the
          Genesets and the column which contains the Genes. If you are unsure
          about the format of the input data, please check out the Welcome
          panel.",
          br(),
          "Upon choosing the respective columns, we will rename the column
          containing the geneset identifiers to 'Genesets' and the column
          containing the gene lists to 'Genes'. ",
          br(),
          selectInput(
            inputId = "alt_name_genesets",
            label = "Which column contains the Geneset identifiers?",
            choices = c(colnames(reactive_values$genesets)),
            multiple = FALSE,
            selected = NULL
          ),
          selectInput(
            inputId = "alt_name_genes",
            label = "Which column contains the Genes?",
            choices = c(colnames(reactive_values$genesets)),
            multiple = FALSE,
            selected = NULL
          ),
          actionButton(
            inputId = "alt_names_start",
            label = "Submit your alternative columns",
            icon = icon("play-circle"),
            style = .actionButtonStyle
          )
        )
      )
    })

    output$dt_genesets <- DT::renderDataTable({
      validate(need(!(
        is.null(reactive_values$genesets)
      ),
      message = "Please provide input data via the button on the left."))

      DT::datatable(
        reactive_values$genesets,
        options = list(scrollX = TRUE, scrollY = "400px"),
        rownames = FALSE
      )
    })

    output$ui_filter_data <- renderUI({
      if (is.null(reactive_values$genesets) ||
          is.null(reactive_values$genes) ||
          is.null(reactive_values$gs_names)) {
        return(NULL)
      }
      tagList(
        box(
          id = "optional_filtering_box",
          width = 12,
          title = "Optional Filtering Step",
          status = "info",
          solidHeader = TRUE,
          h2("Filter your provided Genesets"),
          collapsible = TRUE,
          collapsed = TRUE,
          fluidRow(
            "It might be beneficial to your analysis to filter out general terms
            and genesets before proceeding with the next steps."
          ),
          fluidRow(
            column(
              width = 6,
              sliderInput(
                inputId = "bins_gs_hist",
                label = "Select the bins for the histogram",
                min = 0,
                max = max(vapply(
                  reactive_values$genes, length, numeric(1)
                )),
                value = c(0, max(
                  vapply(reactive_values$genes, length, numeric(1))
                ))
              )
            ),
            column(
              width = 6,
              sliderInput(
                inputId = "bindwidth_hist",
                label = "Select width of the bins",
                min = 1,
                max = 20,
                value = 5
              )
            )
          ),
          fluidRow(column(
            width = 12,
            plotOutput(
              "histogram_initial_data",
              brush = brushOpts(
                "plot_brush",
                resetOnNew = TRUE,
                direction = "x"
              )
            )
          )),
          fluidRow(column(
            width = 12,
            DT::dataTableOutput("table")
          )),
          fluidRow(
            column(
              width = 6,
              selectizeInput(
                inputId = "select_filter_genesets",
                label = "Select individual genesets to be filtered",
                choices = c(reactive_values$gs_description),
                multiple = TRUE
              ),
              actionButton(
                inputId = "filter_genesets",
                label = "Remove the selected Genesets",
                icon = icon("play-circle"),
                style = .actionButtonStyle
              )
            ),
            column(
              width = 6,
              textInput(inputId = "select_filter_genesets_threshold",
                        label = "Filter genesets with size =>"),
              downloadButton(
                outputId = "download_filtered_data",
                label = "Download the filtered data",
                icon = icon("download"),
                style = .actionButtonStyle
              )
            )
          )
        )
      )
    })

    output$histogram_initial_data <- renderPlot({
      gsHistogram(
        reactive_values$genes,
        reactive_values$gs_names,
        gs_description = reactive_values$gs_description,
        start = input$bins_gs_hist[1],
        end = input$bins_gs_hist[2],
        binwidth = input$bindwidth_hist
      )
    })


    output$ui_panel_specify_species <- renderUI({
      if (is.null(reactive_values$genesets) ||
          is.null(reactive_values$genes) ||
          is.null(reactive_values$gs_names)) {
        return(NULL)
      }
      box(
        width = 12,
        title = "Step 2 (optional)",
        status = "warning",
        solidHeader = TRUE,
        tagList(
          h2("Select the species of your data"),
          fluidRow(
            column(width = 6,
                   uiOutput("ui_specify_species")),
            column(
              width = 6,
              "Currently we are only supporting organisms which have a Protein-
              Protein-Interaction (PPI) matrix in the STRING Database. If you
              are not sure if your organisms is represented in this database,
              please follow this link and have a look at the available organisms.
              This can also be used to ensure the correct spelling of your
              organism of interest.",
              br(),
              tags$a(
                href = "https://string-db.org/cgi/download?sessionId=bs0ZXo3WieT0",
                "Click here to get to STRING!",
                target = "_blank"
              )
            )
          )
        )
      )
    })

    output$ui_specify_species <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      fluidRow(column(
        width = 12,
        selectizeInput(
          "species",
          label = "Please select the species of your data.",
          choices = c(
            "",
            "Homo sapiens",
            "Mus musculus",
            "Rattus norvegicus",
            "Arabidopsis thaliana",
            "Saccharomyces cerevisiae",
            "Drosophila melanogaster",
            "Danio rerio",
            "Caenorhabiditis elegans"
          ),
          multiple = FALSE,
          options = list(create = TRUE)
        )
      ))
    })

    output$ui_panel_download_ppi <- renderUI({
      if (input$species == "" ||
          is.na(input$species) || is.null(input$species)) {
        return(NULL)
      }
      # TODO: Handle input that is not species available on STRING
      reactive_values$species <- input$species
      box(
        width = 12,
        title = "Step 3 (optional)",
        status = "success",
        solidHeader = TRUE,
        tagList(
          h2("Download the PPI matrix from STRING"),
          "For more information on the downloaded PPI data,
          please have a look at this panels tour.",
          br(),
          p(),
          fluidRow(column(
            width = 6,
            actionButton(
              "download_ppi",
              label = "Download PPI matrix",
              icon = icon("download"),
              style = .actionButtonStyle
            )
          ),
          column(
            width = 6,
            box(
              id = "PPI_preview",
              width = NULL,
              title = "PPI preview",
              status = "info",
              solidHeader = TRUE,
              collapsible = TRUE,
              collapsed = TRUE,
              fluidRow(column(
                width = 12,
                offset = 0.5,
                DT::dataTableOutput("dt_ppi")
              ))
            )
          )),
          fluidRow(
            column(
              width = 6,
              "Or you can upload an already saved Protein-Protein interaction matrix from your machine.",
              br(),
              p(),
              fileInput(
                inputId = "inputppi",
                label = "Upload a PPI matrix",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/tab-separated-values",
                  "text/plain",
                  ".csv",
                  ".tsv",
                  ".RDS",
                  ".xslx"
                ),
                multiple = FALSE
              )
            )
          ),
          fluidRow(column(width = 12,
                          uiOutput(
                            "ui_save_PPI"
                          )))
        )
      )
    })

    output$dt_ppi <- DT::renderDataTable({
      validate(
        need(!(is.null(
          reactive_values$ppi
        )),
        message = "Please download a Protein-Protein Interaction (PPI) matrix via the button on the left.")
      )
      DT::datatable(
        reactive_values$ppi,
        options = list(scrollX = TRUE, scrollY = "400px"),
        rownames = FALSE
      )
    })

    output$ui_save_PPI <- renderUI({
      if (is.null(reactive_values$ppi)) {
        return(NULL)
      }
      fluidRow(
        column(
          width = 12,
          "In order to avoid lengthy downloads of the PPI matrix in the future,
        you can also now download the PPI matrix and save it to your machine.",
          br(),
          p(),
          downloadButton(
            "save_ppi",
            label = "Save PPI matrix",
            icon = icon("download"),
            style = .actionButtonStyle
          )
        )
      )
    })

    # panel Scores ----------------------------------------------------
    output$ui_panel_scores <- renderUI({
      validate(need(!(
        is.null(reactive_values$genesets)
      ) &
        !(is.null(
          reactive_values$genes
        )) & !(
          is.null(reactive_values$gs_names)
        ),
      message = "Please provide your data first before proceeding."))
      tagList(
        box(
          id = "distance_calc_box",
          width = 12,
          title = "Calculate Distance Scores for your Genesets",
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = FALSE,
          h2("Select the distance score"),
          "Please check out the tour of this panel if you
          want to know more about the individual scores",
          p(),
          fluidRow(
            column(
              width = 6,
              radioButtons(
                inputId = "scoringmethod",
                label = "Select the scoring method for your data",
                choices = c(
                  "pMM",
                  "Kappa",
                  "Jaccard",
                  "Meet-Min",
                  "Sorensen-Dice",
                  "GO Distance"
                ),
                selected = character(0)
              )
            ),
            column(width = 6,
                   uiOutput("ui_score_data"))
          ),
          fluidRow(column(
            width = 12,
            br(),
            p(),
            uiOutput("ui_alpha_parameter")
          ))
        ),
        box(
          id = "distance_scores_box",
          width = 12,
          status = "info",
          title = "Geneset Distance Scores",
          collapsible = TRUE,
          collapsed = FALSE,
          solidHeader = TRUE,
          fluidRow(column(
            width = 12,
            uiOutput("ui_distance_scores_visuals")
          ))
        )
      )
    })


    output$ui_score_data <- renderUI({
      if (is.null(input$scoringmethod) || input$scoringmethod == "") {
        return(NULL)
      }
      fluidRow(
        column(
          width = 12,
          strong("Now you can score your data"),
          br(),
          "Attention: If you have many genesets to score,
           this operation may take some time",
          br(),
          actionButton("score_data",
                       label = "Compute the distances between genesets",
                       style = .actionButtonStyle)
        )
      )
    })

    output$ui_alpha_parameter <- renderUI({
      if (input$scoringmethod == "" ||
          is.null(input$scoringmethod) ||
          input$scoringmethod != "pMM") {
        return(NULL)
      }
      fluidRow(
        column(
          width = 6,
          "Please select how strongly Protein-Protein-Interactions
        should be weighted in the pMM score by setting the scaling factor",
          br(),
          p(),
          sliderInput(
            "alpha",
            label = "Scaling Factor Alpha",
            min = 0,
            max = 1,
            value = 0.5,
            step = 0.01
          )
        )
      )
    })

    output$ui_distance_scores_visuals <- renderUI({
      validate(
        need(length(reactive_values$scores) != 0,
             message = "Please compute the distances between the genesets first in the above box")
      )
      fluidRow(
        column(
          width = 6,
          selectInput(
            "plots_distance_score",
            label = "Choose the Distance Score Results to plot",
            choices = c("", names(reactive_values$scores))
          ),
        ),
        column(
          width = 6,
          br(),
          downloadButton(
            outputId = "download_scores",
            label = "Download the distance scores",
            icon = icon("download"),
            style = .actionButtonStyle
          )

        ),
        fluidRow(column(
          width = 12,
          bs4Dash::tabsetPanel(
            id = "tabsetpanel_scores",
            type = "tabs",
            selected = "Distance Scores Heatmap",
            side = "right",
            tabPanel(title = "Distance Scores Heatmap",
                     withSpinner(
                       plotOutput("scores_heatmap",
                                  height = "800px",
                                  width = "1000px")
                     )),
            tabPanel(title = "Distance Scores Dendrogram",
                     fluidRow(
                       column(
                         width = 2,
                         selectInput(
                           inputId = "cluster_method_dendro",
                           label = "Select a clustering method for the Dendrogram",
                           choices = c("average", "single", "complete",
                                       "median", "centroid"),
                           selected = "average",
                           multiple = FALSE
                         )
                       ),
                       column(width = 10,
                              plotlyOutput("scores_dendro",
                                            height = "800px",
                                            width = "1000px")
                            )
                     )),
            tabPanel(title = "Distance Scores Graph",
                     fluidRow(
                       column(
                         width = 2,
                         sliderInput(
                           "similarityScores",
                           "Distance Threshold:",
                           min = 0,
                           max = 1,
                           value = 0.3,
                           step = 0.05
                         ),
                         selectizeInput(
                           inputId = "scores_graph_search",
                           label = "Search for a specific geneset",
                           choices = c("", reactive_values$gs_description),
                           multiple = TRUE,
                           options = list(
                             create = FALSE,
                             placeholder = "",
                             maxItems = "1"
                           )
                         )
                       ),
                       column(width = 10,
                              withSpinner(
                                visNetworkOutput("scores_Network",
                                                 height = "800px")
                              ))
                     ),
                     fluidRow(
                       box(
                         id = "hub_genes_box",
                         width = 12,
                         title = "Graph metrics",
                         status = "info",
                         solidHeader = TRUE,
                         collapsible = TRUE,
                         collapsed = TRUE,
                         DT::dataTableOutput("dt_graph_metrics")
                       )
                     ))
          )
        ))
      )
    })

    output$scores_heatmap <- renderPlot({
      validate(
        need(length(reactive_values$scores) != 0,
             message = "Please compute the distances between the genesets first in the above box.")
      )
      validate(
        need(input$plots_distance_score != "",
             message = "Please select one of your distance score results to be displayed.")
      )
      distanceHeatmap(reactive_values$scores[[input$plots_distance_score]],
                      chars_limit = 20)
    })


    output$scores_dendro <- renderPlotly({
      validate(
        need(length(reactive_values$scores) != 0,
             message = "Please compute the distances between the genesets first in the above box.")
      )
      validate(
        need(input$plots_distance_score != "",
             message = "Please select one of your distance score results to be displayed.")
      )
      scores <- reactive_values$scores[[input$plots_distance_score]]
      distanceDendro(scores,
                     input$cluster_method_dendro)
    })


    reactive_values$scores_graph <- reactive({
      if (is.null(reactive_values$scores) ||
          length(reactive_values$scores) == 0) {
        showNotification(
          "It seems like you do not have any distance scores. Please score you data first in the Scores panel.",
          type = "warning"
        )
      }
      validate(
        need(length(reactive_values$scores) != 0,
             message = "Please compute the distances between the genesets first in the above box.")
      )
      validate(
        need(input$plots_distance_score != "",
             message = "Please select one of your distance score results to be displayed.")
      )

      scores <- reactive_values$scores[[input$plots_distance_score]]
      adj <- getAdjacencyMatrix(scores,
                                input$similarityScores)
      g <- buildGraph(adj,
                      reactive_values$genesets,
                      reactive_values$gs_description)
      return(g)
    })

    output$scores_Network <- renderVisNetwork({
      if (!any(get.edgelist(reactive_values$scores_graph()) != 0)) {
        showNotification(
          "Please select a larger distance threshold as currently no nodes are connected and the graph cannot be properly rendered.",
          type = "warning"
        )
      } else {
        visNetwork::visIgraph(reactive_values$scores_graph()) %>%
          visNodes(color = list(
            background = "#0092AC",
            highlight = "gold",
            hover = "gold"
          )) %>%
          visEdges(color = list(
            background = "#0092AC",
            highlight = "gold",
            hover = "gold"
          )) %>%
          visOptions(
            highlightNearest = list(
              enabled = TRUE,
              degree = 1,
              hover = TRUE
            ),
            nodesIdSelection = TRUE
          ) %>%
          visExport(name = "distance_scores_network",
                    type = "png",
                    label = "Save Distance Scores graph")
      }
    })

    output$dt_graph_metrics <- DT::renderDataTable({
      dt <- .graphMetricsGenesetsDT(reactive_values$scores_graph(),
                                    reactive_values$genesets)
      DT::datatable(dt,
                    rownames = FALSE,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })

    # panel Graph ------------------------------------------------------
    output$ui_panel_graph <- renderUI({
      validate(need(!(
        is.null(reactive_values$scores) &&
          (length(reactive_values$scores) == 0)
      ),
      message = "Please score you genesets first in the Scores
                    panel."))

      tagList(
        box(
          id = "clustering_selection_box",
          width = 12,
          title = "Select the clustering method",
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = FALSE,
          fluidRow(
            column(
              width = 6,
              radioButtons(
                inputId = "select_clustering",
                label = "Choose one of the following Clustering methods",
                choices = c("Louvain",
                            "Markov",
                            "Fuzzy",
                            "kMeans")
              )
            ),
            column(width = 6,
                   uiOutput("ui_cluster"))
          ),
          fluidRow(column(
            width = 6,
            uiOutput("ui_parameters_clustering")
          ))
        ),
        box(
          id = "cluster_graph_box",
          width = 12,
          status = "info",
          title = "Geneset Cluster Graphs",
          collapsible = TRUE,
          collapsed = FALSE,
          solidHeader = TRUE,
          fluidRow(column(width = 12,
                          uiOutput("ui_cluster_graphs")
                          )
                   )
        ),
        fluidRow(
          bs4Card(
            width = 12,
            id = "card_clustering",
            title = "Clustering graph summaries",
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            closable = FALSE,
            fluidRow(column(
              width = 12,
              DT::dataTableOutput("dt_cluster")
            ))
          )
        )
      )
    })
    
    output$ui_cluster_graphs <- renderUI({
      validate(need(!(
        is.null(reactive_values$cluster) &&
          (length(reactive_values$cluster) == 0)
      ),
      message = "Please cluster your genesets first in the Cluster Graph
                    panel."))
      
      fluidRow(column(
        width = 12,
        bs4Dash::tabsetPanel(
            id = "tabsetpanel_cluster",
            type = "tabs",
            selected = "Geneset Graph",
            side = "right",
            tabPanel(title = "Geneset Graph",
                     fluidRow(
                       column(
                         width = 2,
                         selectInput(
                           inputId = "graphColoring",
                           label = "Color the graph by",
                           choices =
                             if (!(is.null(reactive_values$genesets))) {
                               c(NULL, "Cluster", colnames(dplyr::select_if(
                                 reactive_values$genesets,
                                 is.numeric
                               )))
                             } else {
                               c("Cluster")
                             },
                           multiple = FALSE
                         )
                       ),
                       column(width = 10,
                              withSpinner(
                                visNetworkOutput("cluster_Network",
                                                 height = "700px",
                                                 width = "100%")
                              ))
                     )),
            tabPanel(title = "Cluster-Geneset Bipartite Graph",
                     withSpinner(
                       visNetworkOutput(
                         "cluster_geneset_bipartite_Network",
                         height = "800px",
                         width = "100%"
                       )
                     )),
            tabPanel(title = "Cluster Enrichment Terms Word Cloud",
                     fluidRow(
                       column(
                         width = 2,
                         selectInput(
                           "cluster_nb",
                           "Select a cluster",
                           choices = c(seq_len(length(
                             reactive_values$cluster
                           ))),
                           selected = 1
                         )
                       ),
                       column(
                         width = 12,
                         wordcloud2::wordcloud2Output("enrichment_wordcloud_cluster")
                       )
                     ))
          )
        )
      )
    })

    output$ui_cluster <- renderUI({
      fluidRow(
        column(
          width = 12,
          selectInput(
            "clustering_score_selected",
            label = "Select Distance Scoring Results to use",
            choices = c("", names(reactive_values$scores))
          ),
          br(),
          p(),
          strong("Now you can cluster your data"),
          br(),
          "Attention: If you have many Genesets to cluster,
           this operation may take some time",
          br(),
          actionButton("cluster_data",
                       label = "Cluster the Genesets",
                       style = .actionButtonStyle)
        )
      )
    })

    output$ui_parameters_clustering <- renderUI({
      if (input$select_clustering == "Louvain") {
        uiOutput("ui_louvain")
      } else if (input$select_clustering == "Markov") {
        uiOutput("ui_markov")
      } else if (input$select_clustering == "Fuzzy") {
        uiOutput("ui_fuzzy")
      } else if (input$select_clustering == "kMeans") {
        uiOutput("ui_kMeans")
      }
    })

    output$ui_louvain <- renderUI({
      fluidRow(
        h2(
          "Select the clustering threshold for the Louvain clustering."
        ),
        br(),
        column(
          width = 6,
          sliderInput(
            inputId = "louvain_threshold",
            label = "Select a similarity threshold",
            min = 0,
            max = 1,
            value = 0.3,
            step = 0.05
          )
        )
      )
    })

    output$ui_markov <- renderUI({
      fluidRow(
        h2(
          "Select the clustering threshold for the Markov clustering."
        ),
        br(),
        column(
          width = 6,
          sliderInput(
            inputId = "markov_threshold",
            label = "Select a similarity threshold",
            min = 0,
            max = 1,
            value = 0.3,
            step = 0.05
          )
        )
      )
    })

    output$ui_fuzzy <- renderUI({
      fluidRow(
        h2(
          "Select the clustering thresholds for the Fuzzy clustering."
        ),
        column(
          width = 6,
          sliderInput(
            inputId = "simThreshold",
            label = "Similarity Threshold",
            min = 0,
            max = 1,
            value = 0.3,
            step = 0.1
          ),
          sliderInput(
            inputId = "memThreshold",
            label = "Membership Threshold",
            min = 0,
            max = 1,
            value = 0.5,
            step = 0.1
          ),
          sliderInput(
            inputId = "clustThreshold",
            label = "Clustering Threshold",
            min = 0,
            max = 1,
            value = 0.5,
            step = 0.1
          )
        )
      )
    })

    output$ui_kMeans <- renderUI({
      fluidRow(
        h2("Select the center number for kMeans clustering."),
        br(),
        column(
          width = 6,
          sliderInput(
            inputId = "center",
            label = "Select a number of centers to start with",
            min = 1,
            max = length(reactive_values$gs_names),
            value = 1,
            step = 1
          )
        )
      )
    })


    output$cluster_Network <- renderVisNetwork({
      validate(need(!(is.null(
        reactive_values$cluster
      )),
      message = "Please cluster the genesets first in the above box"))

      graph <- reactive_values$cluster_graph()

      if (!any(get.edgelist(graph) != 0)) {
        showNotification(
          "It seems like you don't have any clusters. Please adapt the similarity Threshold above and re-run the Clustering.",
          type = "warning"
        )
      } else {
        visNetwork::visIgraph(graph) %>%
          visOptions(
            highlightNearest = list(
              enabled = TRUE,
              degree = 1,
              hover = TRUE
            ),
            nodesIdSelection = TRUE,
            selectedBy = list(variable = "cluster", multiple = TRUE)
          ) %>%
          visExport(name = "cluster_network",
                    type = "png",
                    label = "Save Cluster graph")
      }
    })

    reactive_values$cluster_graph <- reactive({
      g <- buildClusterGraph(
        reactive_values$cluster,
        reactive_values$genesets,
        reactive_values$gs_names,
        input$graphColoring,
        reactive_values$gs_description
      )
      return(g)
    })



    output$cluster_geneset_bipartite_Network <- renderVisNetwork({
      validate(need(!(is.null(
        reactive_values$cluster
      )),
      message = "Please cluster you genesets first in the above box"))
      
      graph <- reactive_values$bipartite_graph()

      visNetwork::visIgraph(graph) %>%
        visIgraphLayout(layout = "layout_as_bipartite") %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = TRUE
        ) %>%
        visExport(name = "bipartite_network",
                  type = "png",
                  label = "Save Cluster-Geneset bipartite graph")
    })
    
    
    reactive_values$bipartite_graph <- reactive({
      tryCatch(
        expr = {
          g <- getBipartiteGraph(reactive_values$cluster,
                                 reactive_values$gs_names,
                                 reactive_values$genes)
        },
        error = function(cond) {
          showNotification(
            "It seems like your data does not have any clusters. Please adapt the thresholds and try again.",
            type = "error"
          )
          return(NULL)
        }
      )
      return(g)
    })


    output$dt_cluster <- DT::renderDataTable({
      validate(
        need(
          !(is.null(reactive_values$cluster)),
          message = "It seems like the data has not yet been clustered. Please
        cluster your data first with the box above."
        )
      )
      dt_cluster <- .getClusterDatatable(
        reactive_values$cluster,
        reactive_values$gs_names,
        reactive_values$gs_description
      )
      DT::datatable(dt_cluster,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })

    output$enrichment_wordcloud_cluster <-
      wordcloud2::renderWordcloud2({
        validate(
          need(
            !(is.null(reactive_values$cluster)),
            message = "It seems like the data has not yet been clustered. Please
        cluster your data first with the box above."
          )
        )
        cluster <- as.numeric(input$cluster_nb)
        validate(
          need(
            cluster <= length(reactive_values$cluster),
            message = "It seems like you have selected a number larger than the number of available clusters."
          )
        )
        genesets <- reactive_values$cluster[[cluster]]
        genesets_df <- reactive_values$genesets[genesets, ]

        enrichmentWordcloud(genesets_df)
      })

    # Report panel -----------------------------------------------------------
    output$ui_panel_report <- renderUI({
      tagList(fluidRow(column(width = 11),
                       column(
                         width = 1,
                         actionButton(
                           "tour_bookmarks",
                           label = "",
                           icon = icon("question-circle"),
                           style = .actionButtonStyle
                         )
                       )),
              fluidRow(column(width = 12,
                              uiOutput(
                                "ui_bookmarks"
                              ))),
              fluidRow(bs4Dash::column(
                width = 8,
                offset = 2,
                br(),
                br(),
                tags$a(
                  id = "start_report",
                  class = "btn btn-default shiny-download-link",
                  href = "",
                  target = "_blank",
                  download = NA,
                  icon("download"),
                  "Start the generation of the report"
                )
              )))
    })

    output$ui_bookmarks <- renderUI({
      tagList(fluidRow(
        column(
          width = 6,
          bs4Dash::bs4InfoBoxOutput("infobox_book_genesets",
                                    width = 6),
          h5("Bookmarked genesets"),
          DT::dataTableOutput("bookmarks_genesets"),
          downloadButton(
            "btn_export_genesets",
            label = "",
            class = "biocdlbutton"
          ),
          actionButton(
            "btn_reset_genesets",
            label = "",
            icon = icon("trash"),
            style = .actionButtonStyle
          )
        ),
        column(
          width = 6,
          bs4Dash::bs4InfoBoxOutput("infobox_book_cluster",
                                    width = 6),
          h5("Bookmarked cluster"),
          DT::dataTableOutput("bookmarks_cluster"),
          downloadButton(
            "btn_export_cluster",
            label = "",
            class = "biocdlbutton"
          ),
          actionButton(
            "btn_reset_cluster",
            label = "",
            icon = icon("trash"),
            style = .actionButtonStyle
          )
        )
      ))
    })

    output$infobox_book_genesets <- renderbs4InfoBox({
      bs4Dash::bs4InfoBox(
        title = "Bookmarked genesets",
        value = nrow(reactive_values$bookmarked_genesets),
        icon = icon("bookmark"),
        iconElevation = 3,
        color = "success",
        fill = TRUE,
        width = 12
      )
    })

    output$infobox_book_cluster <- renderbs4InfoBox({
      bs4Dash::bs4InfoBox(
        title = "Bookmarked cluster",
        value = nrow(reactive_values$bookmarked_cluster),
        icon = icon("bookmark"),
        iconElevation = 3,
        color = "info",
        fill = TRUE,
        width = 12
      )
    })

    output$bookmarks_genesets <- DT::renderDataTable({
      DT::datatable(reactive_values$bookmarked_genesets, rownames = FALSE)
    })

    output$bookmarks_cluster <- DT::renderDataTable({
      DT::datatable(reactive_values$bookmarked_cluster, rownames = FALSE)
    })

    output$btn_export_genesets <- downloadHandler(
      filename = function() {
        paste0("GeDi_bookmarked_genesets_",
               gsub(" ", "_", gsub(
                 "-", "", gsub(":", "-", as.character(Sys.time()))
               )),
               ".txt")
      },
      content = function(file) {
        writeLines(text = reactive_values$bookmarked_genesets,
                   con = file)
      }
    )

    output$btn_export_cluster <- downloadHandler(
      filename = function() {
        paste0("GeDi_bookmarked_cluster_",
               gsub(" ", "_", gsub(
                 "-", "", gsub(":", "-", as.character(Sys.time()))
               )),
               ".txt")
      },
      content = function(file) {
        writeLines(text = reactive_values$bookmarked_cluster,
                   con = file)
      }
    )


    output$start_report <- downloadHandler(
      filename = paste0(
        Sys.Date(),
        "_",
        round(stats::runif(1) * 100),
        # for not having all w the same name
        "_GeDiReport.html"
      ),
      content = function(file) {
        # temporarily switch to the temp dir, in case you do not have write permission to the current working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        withProgress(
          rmarkdown::render(
            input = system.file("extdata", "report_GeDi.Rmd", package = "GeDi"),
            output_file = file,
            quiet = TRUE
          ),
          message = "Generating the html report",
          detail = "This can take some time"
        )
      }
    )


    # Observers ---------------------------------------------------------------

    # Data Input Panel --------------------------------------------------------
    observeEvent(input$inputgenesetfile, {
      filename <-
        unlist(strsplit(input$inputgenesetfile$datapath, ".", fixed = TRUE))
      format <- tail(filename, n = 1)

      if (format == "RDS") {
        reactive_values$genesets <-
          readRDS(input$inputgenesetfile$datapath)
      } else if (format == "txt") {
        reactive_values$genesets <- readGenesetsTxt()
      } else if (format == "xlsx") {
        reactive_values$genesets <-
          readxl::read_excel(input$inputgenesetfile$datapath)
      } else {
        showNotification(
          "It seems like your input file has not the right format. Please check the Welcome panel for the right input format and provide your data again.",
          type = "error"
        )
      }

      tryCatch(
        expr = {
          columns <- names(reactive_values$genesets)
          stopifnot(any(columns == "Genesets"))
          reactive_values$gs_names <-
            reactive_values$genesets$Genesets
          reactive_values$gs_description <-
            .getGenesetDescriptions(reactive_values$genesets)
          reactive_values$genes <-
            getGenes(reactive_values$genesets)
          reactive_values$genesets$Genes <-
            vapply(reactive_values$genesets$Genes, function(x)
              gsub("/", ",", x), character(1))
        },
        error = function(cond) {
          showNotification(
            "It seems like your data does not have a columns named 'Genesets' and 'Genes'. Please check your data and try to input it again.",
            type = "error"
          )
          reactive_values$alt_names <- TRUE
          reactive_values$gs_names <- NULL
          reactive_values$genes <- NULL

        }
      )
      # reset all reactive_values in case data has already been loaded before
      reactive_values$ppi <- NULL
      reactive_values$scores <- list()
      reactive_values$seeds <- NULL
      reactive_values$cluster <- NULL
      reactive_values$bookmarked_genesets <- NULL
      reactive_values$bookmarked_cluster <- NULL
    })

    observeEvent(input$alt_names_start, {
      reactive_values$gs_names <-
        reactive_values$genesets[, input$alt_name_genesets]

      reactive_values$genes <- getGenes(reactive_values$genesets,
                                        input$alt_name_genes)

      reactive_values$gs_description <-
        .getGenesetDescriptions(reactive_values$genesets)

      names(reactive_values$genesets)[names(reactive_values$genesets) == input$alt_name_genesets] <-
        "Genesets"
      names(reactive_values$genesets)[names(reactive_values$genesets) == input$alt_name_genes] <-
        "Genes"

      reactive_values$genesets$Genes <-
        vapply(reactive_values$genesets$Genes, function(x)
          gsub("/", ",", x), character(1))

      showNotification(
        "Successfully selected your columns and read your data. We renamed the columns to 'Genesets' and 'Genes'. You can now proceed.",
        type = "message"
      )
      reactive_values$alt_names <- FALSE
    })

    observeEvent(input$btn_loaddemo, {
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Loading demo data now", value = 0)

      data_env <- new.env(parent = emptyenv())
      data("macrophage_topGO_example", envir = data_env, package = "GeDi")
      macrophage_topGO_example <- data_env[["macrophage_topGO_example"]]

      progress$inc(1 / 3, detail = "Extracting Genesets")

      reactive_values$genesets <- macrophage_topGO_example
      reactive_values$gs_names <- macrophage_topGO_example$Genesets
      reactive_values$gs_description <-
        macrophage_topGO_example$Term

      progress$inc(1 / 3, detail = "Extracting Genes")
      reactive_values$genes <- getGenes(reactive_values$genesets)

      progress$inc(1 / 3, detail = "Successfully loaded demo data")
    })

    observeEvent(input$btn_showDataStructure, {
      showModal(
        modalDialog(
          title = "Data structure of the Input Data",
          "This is an example of the required structure for the input data",
          HTML('<img src="GeDi/DataExample.png">'),
          "The data should contain two columns, 'Genesets' and 'Genes', which are tab-seprated.",
          "The list of genes for each genesets should be comma-, tab- or backslash-separated.",
          easyClose = TRUE,
          footer = NULL
        )
      )
    })

    observeEvent(input$plot_brush, {
      df <- buildHistogramData(
        genesets = reactive_values$genes,
        gs_names = reactive_values$gs_names,
        gs_description = reactive_values$gs_description,
        start = input$bins_gs_hist[1],
        end = input$bins_gs_hist[2]
      )

      info_plot <- brushedPoints(df, input$plot_brush, xvar = "Size")
      output$table <- DT::renderDataTable(info_plot)
    })

    observeEvent(input$filter_genesets, {
      if (is.null(input$select_filter_genesets) &&
          input$select_filter_genesets_threshold == "") {
        showNotification("No Genesets selected for filtering.",
                         type = "message")
      } else{
        remove_specific <- list()
        remove_threshold <- list()

        if (!is.null(input$select_filter_genesets)) {
          index <- which(reactive_values$gs_description == input$select_filter_genesets)
          names <- reactive_values$gs_names[index]
          remove_specific <- names
        }

        if (input$select_filter_genesets_threshold != "") {
          threshold <- as.numeric(input$select_filter_genesets_threshold)
          data <- buildHistogramData(reactive_values$genes,
                                     reactive_values$gs_names)
          remove_threshold <- data[data$Size >= threshold,]$Geneset
          if (length(remove_threshold) == 0) {
            showNotification(
              paste(
                "No Genesets with size larger than ",
                threshold,
                " in input data.",
                sep = ""
              ),
              type = "message"
            )
          }
        }

        remove <- c(remove_specific, remove_threshold)

        filtered_data <- .filterGenesets(remove,
                                         reactive_values$genesets)

        reactive_values$filtered_genesets <-
          rbind(reactive_values$filtered_genesets,
                reactive_values$genesets[reactive_values$genesets$Genesets %in% unlist(remove), ])
        reactive_values$genesets <- filtered_data$Geneset
        reactive_values$gs_names <- filtered_data$gs_names
        reactive_values$genes <- filtered_data$Genes
        reactive_values$gs_description <-
          .getGenesetDescriptions(reactive_values$genesets)

        showNotification("Successfully filtered the selected Genesets.",
                         type = "message")
      }
    })

    output$download_filtered_data <- downloadHandler(
      filename = "Filtered_data_GeDi.RDS",
      content = function(file) {
        saveRDS(reactive_values$genesets, file)
      }
    )


    observeEvent(input$download_ppi, {
      validate(need(!(input$species == ""),
                    message = "Please specify the species of your data"))
      if (length(strsplit(input$species, "\\s+")) > 1) {
        showNotification(
          "It seems like you have selected more than one species. Please go back to the box and select the correct species.",
          type = "error"
        )
        return()
      }
      reactive_values$species <- input$species
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Downloading PPI from STRINGdb", value = 0)
      progress$inc(1 / 12, detail = "Get species ID")
      id <- getId(reactive_values$species, cache = TRUE)
      progress$inc(1 / 12, detail = "Validate species ID")
      validate(
        need(
          !is.na(id),
          message = "We could not find your specified species.
                   Please check the spelling and try again."
        )
      )
      progress$inc(1 / 12, detail = "Get species specific STRINGdb")
      stringdb <- getStringDB(as.numeric(id), cache_location = TRUE)
      stringdb

      progress$inc(4 / 12, detail = "Get Annotation information")
      anno_df <- getAnnotation(stringdb)

      progress$inc(1 / 12, detail = "Download PPI")
      reactive_values$ppi <- getPPI(reactive_values$genes,
                                    stringdb = stringdb,
                                    anno_df = anno_df)
      progress$inc(4 / 12, detail = "Successfully downloaded PPI")
    })

    observeEvent(input$inputppi, {
      filename <-
        unlist(strsplit(input$inputppi$datapath, ".", fixed = TRUE))
      format <- tail(filename, n = 1)

      if (format == "RDS") {
        ppi <-
          readRDS(input$inputppi$datapath)
      } else {
        showNotification(
          "It seems like your input file has not the right format. Please check the Welcome panel for the right input format and provide your data again.",
          type = "error"
        )
      }

      reactive_values$ppi <- .checkPPI(ppi)

      # reset all reactive_values in case data has already been loaded before
      reactive_values$scores <- list()
      reactive_values$seeds <- NULL
      reactive_values$cluster <- NULL
      reactive_values$bookmarked_genesets <- NULL
      reactive_values$bookmarked_cluster <- NULL
    })

    output$save_ppi <- downloadHandler(
      filename = function() {
        paste(input$species, "_PPI.RDS", sep = "")
      },
      content = function(file) {
        saveRDS(reactive_values$ppi, file)
      }
    )


    # Scoring panel -----------------------------------------------------------
    observeEvent(input$score_data, {
      if (input$scoringmethod == "" ||
          (is.null(input$scoringmethod))) {
        showNotification(
          "It seems like you did not select a scoring method. Please select a scoring method on the left.",
          type = "error"
        )
      }

      if ((length(reactive_values$genes) == 0 ||
           (is.null(reactive_values$genes)))) {
        showNotification(
          "It seems like the file you've provided does not contain any genesets. Please check you input and retry.",
          type = "error"
        )
      }

      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      reactive_values$cluster <- NULL

      progress$set(message = "Scoring your genesets", value = 0)

      if (input$scoringmethod == "Meet-Min") {
        scores <- getMeetMinMatrix(reactive_values$genes, progress)
      } else if (input$scoringmethod == "Kappa") {
        scores <- getKappaMatrix(reactive_values$genes, progress)
      } else if (input$scoringmethod == "pMM") {
        if (is.null(reactive_values$ppi)) {
          showNotification(
            "It seems like you have not downloaded a PPI matrix. Please return to the Data Input panel and download the respective PPI.",
            type = "error"
          )
          scores <- NULL
        } else {
          scores <- getpMMMatrix(
            reactive_values$genes,
            reactive_values$ppi,
            alpha = input$alpha,
            progress = progress
          )
        }
      } else if (input$scoringmethod == "Jaccard") {
        scores <- getJaccardMatrix(reactive_values$genes,
                                   progress = progress)
      } else if (input$scoringmethod == "Sorensen-Dice") {
        scores <- getSorensenDiceMatrix(reactive_values$genes,
                                        progress = progress)
      } else if (input$scoringmethod == "GO Distance") {
        tryCatch(
          expr = {
            scores <- goDistance(reactive_values$gs_names,
                                   progress = progress)
          },
          error = function(cond) {
            showNotification(
              paste(
                "It seems like there occured the following error: ",
                cond,
                paste = ""
              ),
              type = "error"
            )
            scores <- NULL
          }
        )
      }

      if (is.null(scores)) {
        showNotification(
          "It seems like something went wrong while scoring your data. Most likely this is due to the genesets not containg genes.",
          type = "error"
        )
      } else {
        rownames(scores) <- colnames(scores) <- reactive_values$gs_names
        l <- length(reactive_values$scores) + 1

        reactive_values$scores[[l]] <- scores
        names(reactive_values$scores)[[l]] <- input$scoringmethod
        updateBox("distance_calc_box", action = "toggle")
      }
    })

    output$download_scores <- downloadHandler(
      filename = function() {
        paste(input$scoringmethod, "_GeDi.RDS", sep = "")
      },
      content = function(file) {
        saveRDS(reactive_values$scores[[input$scoringmethod]], file)
      }
    )

    observeEvent(input$scores_graph_search, {
      current_node <- input$scores_graph_search
      if (current_node != "") {
        isolate({
          if (current_node %in% reactive_values$gs_description) {
            current_node <-
              reactive_values$gs_names[which(reactive_values$gs_description == current_node)]
            visNetworkProxy("scores_Network") %>% visSelectNodes(id = current_node)
          }
        })
      } else {
        showNotification(
          "It seems like the geneset identifier you searched for cannot be found in your input data. Please select a different geneset.",
          type = "message"
        )
      }
    })


    # Clustering panel --------------------------------------------------------
    observeEvent(input$cluster_data, {
      if (is.null(reactive_values$scores) ||
          length(reactive_values$scores) == 0) {
        showNotification(
          "It seems like you did not compute the distances between the genesets yet. Please go back to the Distance Scores panel and select a score of your choice.",
          type = "error"
        )
      }else if (input$clustering_score_selected == "" || length(input$clustering_score_selected) == 0) {
        showNotification(
          "Please select a score of your choice to start the clustering.",
          type = "error"
        )
      } else {
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      scores <-
        reactive_values$scores[[input$clustering_score_selected]]

      if (input$select_clustering == "Louvain") {
        progress$set(message = "Start the Louvain Clustering", value = 0)
        cluster <- clustering(scores,
                              input$louvain_threshold,
                              cluster_method = "louvain")
        progress$inc(0.8, detail = "Finished Louvain clustering")
      } else if (input$select_clustering == "Markov") {
        progress$set(message = "Start the Markov Clustering", value = 0)
        cluster <- clustering(scores,
                              input$markov_threshold,
                              cluster_method = "markov")
        progress$inc(0.8, detail = "Finished Markov clustering")
      } else if (input$select_clustering == "Fuzzy") {
        progress$set(message = "Start the fuzzy clustering", value = 0)
        progress$inc(0.1, detail = "Finding intial seeds")

        seeds <- seedFinding(scores,
                             input$simThreshold,
                             input$memThreshold)

        progress$inc(0.4,
                     detail = "Found initial seeds. Now clustering the data.")

        cluster <- fuzzyClustering(seeds, input$clustThreshold)

        progress$inc(0.3, detail = "Finished clustering the data.")

        if (is.null(seeds)) {
          showNotification(
            "It seems like something went wrong while clustering your data. Most likely this is due to the genesets not containg genes.",
            type = "error"
          )
        } else{
          reactive_values$seeds <- seeds
        }
      } else if (input$select_clustering == "kMeans") {
        progress$set(message = "Start the kMeans Clustering", value = 0)
        cluster <- kMeansClustering(scores,
                                  input$center)
        progress$inc(0.8, detail = "Finished clustering the data.")
      }

      if (is.null(cluster)) {
        showNotification(
          "It seems like something went wrong while clustering your data. Most likely this is due to an error in the scoring.",
          type = "error"
        )
      } else{
        reactive_values$cluster <- cluster
        progress$inc(0.2, detail = "Successfully clustered data.")
        updateBox("clustering_selection_box", action = "toggle")
      }
    }})



    # Report panel ----------------------------------------------------------

    observeEvent(input$bookmarker, {
      if (input$tabs == "tab_welcome" ||
          input$tabs == "tab_data_input" ||
          input$tabs == "tab_scores") {
        showNotification("Welcome to GeDi! There is nothing to Bookmark here.")
      } else if (input$tabs == "tab_graph") {
        if (input$tabsetpanel_cluster == "Geneset Graph") {
          g <- reactive_values$cluster_graph()
          cur_sel <- input$cluster_Network_selected
          if (cur_sel == "") {
            showNotification("Select a node in the network to bookmark it",
                             type = "warning")
          } else{
            cur_node <- match(cur_sel, V(g)$name)
            cur_sel_term <- reactive_values$gs_description[cur_node]
            cur_sel_merged <-
              c("Geneset_id" = cur_sel,
                "Geneset_description" = cur_sel_term)
            if (cur_sel %in% reactive_values$bookmarked_genesets) {
              showNotification(
                sprintf(
                  "The selected gene set %s is already in the set of the bookmarked genesets.",
                  cur_sel
                ),
                type = "default"
              )
            } else {
              reactive_values$bookmarked_genesets <-
                unique(rbind(
                  reactive_values$bookmarked_genesets,
                  cur_sel_merged
                ))
              showNotification(
                sprintf(
                  "Added %s to the bookmarked genesets. The list contains now %d elements",
                  cur_sel,
                  nrow(reactive_values$bookmarked_genesets)
                ),
                type = "message"
              )
            }
          }
        }
        else if (input$tabsetpanel_cluster == "Cluster-Geneset Bipartite Graph") {
          g <- reactive_values$bipartite_graph()
          cur_sel <-
            input$cluster_geneset_bipartite_Network_selected
          cur_node <- match(cur_sel, V(g)$name)
          cur_nodetype <- V(g)$nodeType[cur_node]

          if (cur_sel == "") {
            showNotification("Select a node in the network to bookmark it",
                             type = "warning")
          } else {
            if (cur_nodetype == "Geneset") {
              cur_sel_term <- reactive_values$gs_description[cur_node]
              cur_sel_merged <-
                c("Geneset_id" = cur_sel,
                  "Geneset_description" = cur_sel_term)
              if (cur_sel %in% reactive_values$bookmarked_genesets) {
                showNotification(
                  sprintf(
                    "The selected gene set %s is already in the set of the bookmarked genesets.",
                    cur_sel
                  ),
                  type = "default"
                )
              } else {
                reactive_values$bookmarked_genesets <-
                  unique(rbind(
                    reactive_values$bookmarked_genesets,
                    cur_sel_merged
                  ))
                showNotification(
                  sprintf(
                    "Added %s to the bookmarked genesets. The list contains now %d elements",
                    cur_sel,
                    nrow(reactive_values$bookmarked_genesets)
                  ),
                  type = "message"
                )
              }
            } else if (cur_nodetype == "Cluster") {
              cluster_nb <- as.numeric(unlist(strsplit(cur_sel, " "))[[2]])
              cur_sel_cluster <-
                reactive_values$cluster[[cluster_nb]]
              cur_sel_cluster_member <-
                reactive_values$gs_names[cur_sel_cluster]
              cur_sel_cluster_term <-
                reactive_values$gs_description[cur_sel_cluster]
              cur_cluster <- c(
                "Cluster" = cur_sel,
                "Cluster_Members" = paste(cur_sel_cluster_member,
                                          collapse = ", "),
                "Cluster_Member_Description" = paste(cur_sel_cluster_term,
                                                     collapse = ", ")
              )
              if (cur_sel %in% reactive_values$bookmarked_cluster) {
                showNotification(
                  sprintf(
                    "The selected cluster %sis already in the set of the bookmarked cluster.",
                    cur_sel
                  ),
                  type = "default"
                )
              } else {
                reactive_values$bookmarked_cluster <-
                  unique(rbind(
                    reactive_values$bookmarked_cluster,
                    cur_cluster
                  ))
                showNotification(
                  sprintf(
                    "Added %s to the bookmarked cluster. The list contains now %d elements",
                    cur_sel,
                    nrow(reactive_values$bookmarked_cluster)
                  ),
                  type = "message"
                )
              }
            }
          }
        }
      }
      else if (input$tabs == "tab_report") {
        showNotification(
          "You are already in the Report tab where you can see your bookmarked genesets and cluster."
        )
      }
    })

    observeEvent(input$btn_reset_genesets, {
      reactive_values$bookmarked_genesets <- c()
      showNotification("Resetting list of bookmarked genesets...")
    })

    observeEvent(input$btn_reset_cluster, {
      reactive_values$bookmarked_cluster <- c()
      showNotification("Resetting list of bookmarked cluster...")
    })

    # Tour Observers ----------------------------------------------------------
    observeEvent(input$tour_welcome, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_welcome.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_data_input, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_data_input.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_scoring, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_scoring.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session,
              options = list(steps = tour, tooltipPosition = "bottom"))
    })

    observeEvent(input$tour_graph, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_graph.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_bookmarks, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_bookmarks.txt",
                    package = "GeDi"),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$btn_info_session, {
      showModal(
        modalDialog(
          title = "Session information",
          size = "l",
          fade = TRUE,
          footer = NULL,
          easyClose = TRUE,
          tagList(tags$code("> sessionInfo()"),
                  renderPrint({
                    utils::sessionInfo()
                  }))
        )
      )
    })

    observeEvent(input$btn_info_gedi, {
      showModal(
        modalDialog(
          title = "About GeDi",
          size = "l",
          fade = TRUE,
          footer = NULL,
          easyClose = TRUE,
          tagList(includeMarkdown(
            system.file("extdata", "about.md", package = "GeDi")
          ),
          renderPrint({
            utils::citation("GeDi")
          }))
        )
      )
    })
  }
  # nocov end

  shinyApp(ui = gedi_ui, server = gedi_server)
}

