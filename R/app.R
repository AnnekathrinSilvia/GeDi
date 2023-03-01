#' GeDi main function
#'
#' TODO: In the end add to all the examples
#' also some with the example data not only with some to see the data structure
#'
#' @param genesets a dataframe of genesets containing the columns: genesets (the names / ids of the sets) and genes (the genes included in the genesets)
#' @param ppi a Protein-Protein interaction matrix
#' @param alpha a scaling factor in between 0 and 1
#'
#' @return A Shiny app object is returned
#' @export
#' @import visNetwork
#' @import shiny
#' @import plotly
#' @import shinyBS
#' @import fontawesome
#' @import bs4Dash
#' @importFrom rintrojs introjs
#' @importFrom utils read.delim data
#' @importFrom shinycssloaders withSpinner
#' @importFrom igraph V degree delete_vertices
#'
#'
#' @examples
#' \dontrun{
#' GeDi()
#'
#' # Alternatively, you can also start the application with your data directly
#' loaded.
#'
#' data(macrophage_topGO_example, package = "GeDi")
#' GeDi(genesets = macrophage_topGO_example)
#' }
#'
GeDi <- function(genesets = NULL,
                 ppi = NULL,
                 alpha = 1) {
  oopt <- options(spinner.type = 6, spinner.color = "#0092AC")
  on.exit(options(oopt))

  usage_mode <- "shiny_mode"

  # UI definition -----------------------------------------------------------

  # dashpage definition -----------------------------------------------------
  gedi_ui <- bs4DashPage(
    title = "GeDi",
    dark = NULL,
    # navbar definition -------------------------------------------------------
    header = bs4DashNavbar(
      tagList(tags$code(tags$h3("GeDi"))),
      tags$span(style = "display:inline-block; width: 30%"),
      tagList(
        shinyWidgets::dropdownButton(
          inputId = "ddbtn_docs",
          circle = FALSE,
          status = "info",
          icon = icon("book"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "More documentation on GeDi"),
          tags$h5("Documentation"),
          actionButton(
            inputId = "btn_first_help",
            icon = icon("circle-question"),
            label = "First Help",
            style = .actionButtonStyle
          )
        ),
        shinyWidgets::dropdownButton(
          inputId = "ddbtn_info",
          circle = FALSE,
          status = "info",
          icon = icon("info"),
          width = "300px",
          size = "xs",
          right = TRUE,
          tooltip = shinyWidgets::tooltipOptions(title = "More info on GeDi and on the current session"),
          tags$h5("Additional information"),
          actionButton(
            inputId = "btn_docs_vignette",
            icon = icon("book-open"),
            label = "Open GeDi Vignette", style = .actionButtonStyle,
            onclick = ifelse(system.file("doc", "GeDi_manual.html", package = "GeDi") != "",
                             "",
                             "window.open('https://federicomarini.github.io/GeneTonic/articles/GeneTonic_manual.html', '_blank')"
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
      )
    ),

    # sidebar definition ------------------------------------------------------
    sidebar = bs4DashSidebar(
      id = "sidebar",
      title = HTML("<small>GeDi</small>"),
      # src = "GeneTonic/GeneTonic.png",
      skin = "dark",
      status = "info",
      brandColor = NULL,
      # url = "https://bioconductor.org/packages/GeneTonic",
      collapsed = TRUE,
      elevation = 1,
      opacity = 0.8,
      bs4SidebarMenu(
        bs4SidebarMenuItem("Welcome!",
          tabName = "tab_welcome",
          icon = icon("house")
        ),
        bs4SidebarMenuItem(
          "Data Upload",
          tabName = "tab_data_upload",
          icon = icon("square-share-nodes")
        ),
        bs4SidebarMenuItem(
          "Distance Scores",
          tabName = "tab_scores",
          icon = icon("diagram-project")
        ),
        bs4SidebarMenuItem(
          "Graph",
          tabName = "tab_graph",
          icon = icon("square-share-nodes")
        ),
        bs4SidebarMenuItem(
          "Report",
          tabName = "tab_report",
          icon = icon("file")
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
        bs4TabItem(
          tabName = "tab_welcome",
          tagList(
            fluidRow(
              column(width = 11),
              column(
                width = 1,
                actionButton(
                  "tour_welcome",
                  label = "",
                  icon = icon("circle-question"),
                  style = .tourButtonStyle
                )
              )
            ),
            uiOutput("ui_panel_welcome")
          )
        ),
        # ui panel data upload -----------------------------------------------
        bs4TabItem(
          tabName = "tab_data_upload",
          tagList(
            fluidRow(
              column(width = 11),
              column(
                width = 1,
                br(),
                actionButton(
                  "tour_data_upload",
                  label = "",
                  icon = icon("circle-question"),
                  style = .tourButtonStyle
                ),
                shinyBS::bsTooltip(
                  id = "tour_data_upload",
                  title = "Click me to start a tour of this section!",
                  placement = "top",
                  trigger = "hover",
                  options = list(container = "body")
                ),
                style = "float:right"
              )
            ),
            uiOutput("ui_panel_data_upload"),
            uiOutput("ui_filter_data"),
            uiOutput("ui_panel_specify_species"),
            uiOutput("ui_panel_download_ppi")
          )
        ),
        # ui panel scores -------------------------------------------------
        bs4TabItem(
          tabName = "tab_scores",
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
                  style = .tourButtonStyle
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
          )
        ),
        # ui panel graph ---------------------------------------------------
        bs4TabItem(
          tabName = "tab_graph",
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
                  style = .tourButtonStyle
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
          )
        ),
        bs4TabItem(
          tabName = "tab_report",
          uiOutput("ui_panel_report")
        )
      )
    ),
    # controlbar definition ------------------------------------------------
    controlbar = bs4Dash::bs4DashControlbar(
      collapsed = TRUE,
      icon = icon("gears"),
      uiOutput("ui_controlbar")
    ),

    # footer definition -------------------------------------------------------
    footer = bs4Dash::bs4DashFooter(
      left =
        fluidRow(
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
            "- Institute for Medical Biostatistics, Epidemiology and
            Informatics",
            br(),
            "License: ",
            tags$a(href = "https://opensource.org/licenses/MIT", "MIT"),
            "- The GeDi package is developed and available on ",
            tags$a(href = "https://github.com/AnnekathrinSilvia/GeDi", "GitHub")
          )
        ),
      right = NULL
    )
  )
  # server deifnition ---------------------------------------------------------
  gedi_server <- function(input, output, session) {
    # initializing reactives --------------------------------------------------
    reactive_values <- reactiveValues()
    reactive_values$genesets <- NULL
    reactive_values$gs_names <- NULL
    reactive_values$genes <- NULL
    reactive_values$species <- NULL
    reactive_values$ppi <- NULL
    reactive_values$scores <- NULL
    reactive_values$seeds <- NULL
    reactive_values$cluster <- NULL
    reactive_values$alt_names <- FALSE


    # panel Welcome ----------------------------------------------------------
    output$ui_panel_welcome <- renderUI({
      tagList(fluidRow(
        column(
          width = 12,
          includeMarkdown(system.file("extdata", "welcome.md", package = "GeDi")),
          br(), br(),
          p("If you see a grey box like this one open below..."),

          shinyBS::bsCollapse(
            id = "help_welcome", open = "Help",
            shinyBS::bsCollapsePanel(
              "Help",
              includeMarkdown(system.file("extdata", "help_welcome.md", package = "GeDi"))
            )
          ),

          actionButton("introexample", "If you see a button like this...", icon("info"),
                       style = "color: #ffffff; background-color: #0092AC; border-color: #2e6da4"
          ),
          p("... you can click on that to start a tour based on introJS"),
          br(), br()
        )
      )
      )
    })


    # panel Data Upload ------------------------------------------------------
    output$ui_panel_data_upload <- renderUI({
      tagList(
        box(
          width = 12,
          title = "Step 1",
          status = "danger",
          solidHeader = TRUE,
          h2("Upload your Genesets input data"),
          fluidRow(
            column(
              width = 6,
              fileInput(
                inputId = "uploadgenesetfile",
                label = "Upload a geneset file",
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
              ),
              br(),
              "... or you can also ",
              actionButton(
                "btn_loaddemo",
                "Load the demo data",
                icon = icon("play-circle"),
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
                id = "Genests_preview",
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
      if (is.null(input$uploadgenesetfile)) {
        return(NULL)
      }

      guessed_sep <-
        sepguesser(input$uploadgenesetfile$datapath)
      genesets <-
        utils::read.delim(
          input$uploadgenesetfile$datapath,
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
          "It seems like your data does not contain columns named Genesets and Genes.
        Please select the column which contains the identifiers for the Genesets
        and the column which contains the Genes. If you are unsure about the
        format of the input data, please check out the Welcome panel.",
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
      validate(need(
        !(
          is.null(reactive_values$genesets)
        ),
        message = "Please upload a text file via the button on the left."
      ))

      DT::datatable(reactive_values$genesets,
        options = list(scrollX = TRUE, scrollY = "400px")
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
          h2("Filter your uploaded Genesets"),
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
                max = max(sapply(reactive_values$genes, length)),
                value = c(0, max(
                  sapply(reactive_values$genes, length)
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
                resetOnNew = T,
                direction = "x"
              )
            )
          )),
          fluidRow(column(
            width = 12,
            DT::dataTableOutput("table")
          )),
          fluidRow(column(
            width = 6,
            selectizeInput(
              inputId = "select_filter_genesets",
              label = "Select the genesets to be filtered",
              choices = c(reactive_values$gs_names),
              multiple = TRUE,
              options = list(create = TRUE)
            ),
            actionButton(
              inputId = "filter_genesets",
              label = "Remove the selected Genesets",
              icon = icon("play-circle"),
              style = .actionButtonStyle
            )
          ))
        )
      )
    })

    output$histogram_initial_data <- renderPlot({
      gs_histogram(
        reactive_values$genes,
        reactive_values$gs_names,
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
            column(
              width = 6,
              uiOutput("ui_specify_species")
            ),
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
                href="https://string-db.org/cgi/download?sessionId=bs0ZXo3WieT0",
                "Click here to get to STRING!",
                target = "_blank")
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
            "Homo Sapiens",
            "Mus musculus",
            "Rattus norvegicus",
            "Arabidopsis thaliana",
            "Saccharomyces cerevisiae",
            "Drosophila melanogaster",
            "Danio rerio",
            "Caenorhabiditis elegans"
          ),
          multiple = TRUE,
          options = list(create = TRUE)
        )
      ))
    })

    output$ui_panel_download_ppi <- renderUI({
      if (input$species == "" ||
        is.na(input$species) || is.null(input$species)) {
        return(NULL)
      }
      if (length(strsplit(input$species, "\\s+")) > 1) {
        showNotification(
          "It seems like you have selected more than one species. Please go back to the box and select the correct species.",
          type = "error"
        )
      }
      reactive_values$species <- input$species
      box(
        width = 12,
        title = "Step 3 (optional)",
        status = "success",
        solidHeader = TRUE,
        tagList(
          h2("Download the PPI matrix from STRING"),
          fluidRow(
            column(
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
            )
          )
        )
      )
    })

    output$dt_ppi <- DT::renderDataTable({
      validate(
        need(
          !(is.null(
            reactive_values$ppi
          )),
          message = "Please download a Protein-Protein Interaction (PPI) matrix via the button on the left."
        )
      )
      DT::datatable(reactive_values$ppi,
        options = list(scrollX = TRUE, scrollY = "400px")
      )
    })

    # panel Scores ----------------------------------------------------
    output$ui_panel_scores <- renderUI({
      tagList(
        box(
          id = "distance_calc_box",
          width = 12,
          title = "Calculate distance scores for you Genesets",
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = FALSE,
          h2("Select the distance score"),
          fluidRow(
            column(
              width = 6,
              radioButtons(
                inputId = "scoringmethod",
                label = "Select the scoring method for your data",
                choices = c("PMM", "Kappa", "Jaccard", "Meet-Min"),
                selected = character(0)
              )
            ),
            column(
              width = 6,
              uiOutput("ui_score_data")
            )
          ),
          fluidRow(column(width = 12))
        ),
        box(
          id = "distance_scores_box",
          width = 12,
          status = "info",
          title = "Geneset Distance Scores",
          collapsible = TRUE,
          collapsed = FALSE,
          solidHeader = TRUE,
          fluidRow(
            column(
              width = 12,
              uiOutput("ui_distance_scores_visuals")
            )
          )
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
            label = "Score the Genesets",
            style = .actionButtonStyle
          )
        )
      )
    })

    output$ui_distance_scores_visuals <- renderUI({
      validate(need(
        !(is.null(
          reactive_values$scores
        )),
        message = "Please score your genesets first in the above box"
      ))
      fluidRow(
        column(
          width = 12,
            bs4Dash::tabsetPanel(
              id = "tabsetpanel_scores",
              type = "tabs",
              selected = "Distance Scores Heatmap",
              side = "right",
              tabPanel(
                title = "Distance Scores Heatmap",
                actionButton(inputId = "create_heatmap",
                             label = "Calculate Distance Score Heatmap",
                             icon = icon("gears"),
                             style = .actionButtonStyle),
                htmlOutput("scores_heatmap")
              ),
              tabPanel(
                title = "Distance Scores Dendrogram",
                withSpinner(
                  plotlyOutput("scores_dendro",
                               height = "800px",
                               width = "1000px"
                  )
                )
              ),
              tabPanel(
                title = "Distance Scores Graph",
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
                      choices = c("", reactive_values$gs_names),
                      multiple = TRUE,
                      options = list(
                        create = FALSE,
                        placeholder = "",
                        maxItems = "1"
                      )
                    )
                  ),
                  column(
                    width = 9,
                    withSpinner(
                      visNetworkOutput("scores_Network",
                                       width = "800px",
                                       height = "800px"
                      )
                    )
                  )
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
                    bs4Dash::tabsetPanel(
                      id = "tabsetpanel_graph_metrics",
                      type = "tabs",
                      selected = "Node degree",
                      side = "right",
                      tabPanel(
                        title = "Node degree",
                        DT::dataTableOutput("dt_degree")
                      ),
                      tabPanel(
                        title = "Node betweenness",
                        DT::dataTableOutput("dt_betweenness")
                      ),
                      tabPanel(
                        title = "Harmonic centrality",
                        DT::dataTableOutput("dt_harmonic_centrality")
                      ),
                      tabPanel(
                        title = "Clustering Coefficient",
                        DT::dataTableOutput("dt_clustering_coefficient")
                      )
                    )
                  )
                )
              )
            )
        )
      )
    })

    scores_heatmap_react <- reactive({
      validate(need(
        !(is.null(
          reactive_values$scores
        )),
        message = "Please score your genesets first in the above box"
      ))
      res <- ComplexHeatmap::draw(distance_heatmap(reactive_values$scores,
                                                   chars_limit = 20,
                                                   hcluster = input$hcluster
      ))
      return(res)
    })


    output$scores_dendro <- renderPlotly({
      distance_dendro(
        reactive_values$scores,
        input$cluster_method_dendro
      )
    })


    reactive_values$scores_graph <- reactive({
      # TODO: Handle empty scores matrix (matrix of nrow() == 0)
      adj <- getAdjacencyMatrix(
        reactive_values$scores,
        input$similarityScores
      )
      g <- buildGraph(adj)
      return(g)
    })

    output$scores_Network <- renderVisNetwork({
      if (!any(igraph::get.edgelist(reactive_values$scores_graph()) != 0)) {
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
          visExport(
            name = "distance_scores_network",
            type = "png",
            label = "Save Distance Scores graph"
          )
      }
    })

    output$dt_degree <- DT::renderDataTable({
      degree <- igraph::degree(reactive_values$scores_graph(),
                               mode = "all")
      dt <- .graphMetricsGenesetsDT(reactive_values$scores_graph(),
                                    reactive_values$genesets,
                                    degree,
                                    "Degree")
      DT::datatable(dt,
                    rownames = FALSE,
                    options = list(scrollX = TRUE, scrollY = "400px")
      )
    })

    output$dt_betweenness <- DT::renderDataTable({
      validate(
        need(
          !(is.null(
            reactive_values$scores
          )),
          message = "Please score you genesets first in the above box"
        )
      )
      betweenness <- igraph::betweenness(reactive_values$scores_graph(),
                                         directed = FALSE)
      dt <- .graphMetricsGenesetsDT(reactive_values$scores_graph(),
                                    reactive_values$genesets,
                                    betweenness,
                                    "Betweenness")
      DT::datatable(dt,
                    rownames = FALSE,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })

    output$dt_harmonic_centrality <- DT::renderDataTable({
      validate(
        need(
          !(is.null(
            reactive_values$scores
          )),
          message = "Please score you genesets first in the above box"
        )
      )
      centrality <- igraph::harmonic_centrality(reactive_values$scores_graph(),
                                                mode = "all")
      dt <- .graphMetricsGenesetsDT(reactive_values$scores_graph(),
                                    reactive_values$genesets,
                                    centrality,
                                    "Harmonic Centrality")
      DT::datatable(dt,
                    rownames = FALSE,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })

    output$dt_clustering_coefficient <- DT::renderDataTable({
      validate(
        need(
          !(is.null(
            reactive_values$scores
          )),
          message = "Please score you genesets first in the above box"
        )
      )
      clustering_coef <- igraph::transitivity(reactive_values$scores_graph(),
                                              type = "global")
      dt <- .graphMetricsGenesetsDT(reactive_values$scores_graph(),
                                    reactive_values$genesets,
                                    clustering_coef,
                                    "Clustering Coefficient")
      DT::datatable(dt,
                    rownames = FALSE,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })


    # panel Graph ------------------------------------------------------

    output$ui_panel_graph <- renderUI({
      validate(need(!is.null(reactive_values$scores),
        message = "Please score you genesets first in the Scores panel."
      ))

      tagList(
        box(
          id = "clustering_param_box",
          width = 12,
          title = "Set the parameters for Clustering",
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = FALSE,
          fluidRow(
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
            ),
            column(
              width = 6,
              uiOutput("ui_cluster")
            )
          ),
          fluidRow(column(width = 12))
        ),
        fluidRow(column(
          width = 12,
          bs4Dash::bs4Card(
            id = "tabcard_cluster",
            title = "Geneset Cluster Graph",
            elevation = 1,
            width = 12,
            closable = TRUE,
            bs4Dash::tabsetPanel(
              id = "tabsetpanel_cluster",
              type = "tabs",
              selected = "Geneset Graph",
              side = "right",
              tabPanel(
                title = "Geneset Graph",
                withSpinner(
                  visNetworkOutput("cluster_Network",
                    height = "700px",
                    width = "100%"
                  )
                )
              ),
              tabPanel(
                title = "Cluster-Geneset Bipartite Graph",
                withSpinner(
                  visNetworkOutput(
                    "cluster_geneset_bipartite_Network",
                    height = "700px",
                    width = "100%"
                  )
                )
              )
            )
          )
        )),
        fluidRow(
          bs4Dash::bs4Card(
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

    output$ui_cluster <- renderUI({
      fluidRow(
        column(
          width = 12,
          strong("Now you can cluster your data"),
          br(),
          "Attention: If you have many Genesets to cluster,
           this operation may take some time",
          br(),
          actionButton("cluster_data",
            label = "Cluster the Genesets",
            style = .actionButtonStyle
          )
        )
      )
    })

    output$cluster_Network <- renderVisNetwork({
      validate(need(
        !(is.null(
          reactive_values$cluster
        )),
        message = "Please cluster you genesets first in the above box"
      ))


      if (!any(igraph::get.edgelist(reactive_values$cluster_graph()) != 0)) {
        showNotification(
          "It seems like you don't have any clusters. Please adapt the similarity Threshold above and re-run the Clustering.",
          type = "warning"
        )
      } else {
        visNetwork::visIgraph(reactive_values$cluster_graph()) %>%
          # visNodes(color = list(background = "#0092AC", highlight = "gold", hover = "gold")) %>%
          # visEdges(color = list(background = "#0092AC", highlight = "gold", hover = "gold")) %>%
          visOptions(
            highlightNearest = list(
              enabled = TRUE,
              degree = 1,
              hover = TRUE
            ),
            nodesIdSelection = TRUE
          ) %>%
          visExport(
            name = "cluster_network",
            type = "png",
            label = "Save Cluster graph"
          )
      }
    })

    reactive_values$cluster_graph <- reactive({
      g <- buildClusterGraph(
        reactive_values$cluster,
        reactive_values$genesets,
        reactive_values$gs_names,
        input$graphColoring
      )
      return(g)
    })



    output$cluster_geneset_bipartite_Network <- renderVisNetwork({
      validate(need(
        !(is.null(
          reactive_values$cluster
        )),
        message = "Please cluster you genesets first in the above box"
      ))

      visNetwork::visIgraph(reactive_values$bipartite_graph()) %>%
        visIgraphLayout(layout = "layout_as_bipartite") %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = TRUE
        ) %>%
        visExport(
          name = "bipartite_network",
          type = "png",
          label = "Save Cluster-Geneset bipartite graph"
        )
    })

    reactive_values$bipartite_graph <- reactive({
      # TODO: Handle no clusters
      g <- getBipartiteGraph(
        reactive_values$cluster,
        reactive_values$gs_names,
        reactive_values$genes
      )
      # g <- add_layout_(as_bipartite(g))
      return(g)
    })


    output$dt_cluster <- DT::renderDataTable({
      validate(
        need(
          !(is.null(
            reactive_values$cluster
          )),
          message = "It seems like the data has not yet been clustered. Please cluster your data first with the box above."
        )
      )
      dt_cluster <- getClusterDatatable(
        reactive_values$cluster,
        reactive_values$gs_names
      )
      DT::datatable(dt_cluster,
        options = list(scrollX = TRUE, scrollY = "400px")
      )
    })

    # controlbar --------------------------------------------------------------

    output$ui_controlbar <- renderUI({
      tagList(
        checkboxInput("hcluster",
          label = "Activate hclust on the heatmap",
          value = FALSE
        ),
        selectInput(
          inputId = "cluster_method_dendro",
          label = "Select a clustering method for the Dendrogram",
          choices = c(
            "average", "single", "complete",
            "median", "centroid"
          ),
          selected = "average",
          multiple = FALSE
        ),
        numericInput(
          inputId = "n_genesets",
          label = "Number of genesets",
          value = 15,
          min = 1,
          max = 100
        ),
        selectInput(
          "graphColoring",
          label = "Color the nodes by: ",
          choices = if (!(is.null(reactive_values$genesets))) {
            c(NULL, colnames(
              dplyr::select_if(reactive_values$genesets, is.numeric)
            ))
          } else {
            c(NULL)
          },
          selected = NULL,
          multiple = FALSE
        )
      )
    })
    outputOptions(output, "ui_controlbar", suspendWhenHidden = FALSE)


    # Report panel -------------------------------------------------------------
    output$ui_panel_report <- renderUI({
      tagList(
        fluidRow(
          column(
            width = 11
          ),
          column(
            width = 1,
            actionButton(
              "tour_bookmarks",
              label = "", icon = icon("question-circle"),
              style = .tourButtonStyle
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            uiOutput("ui_bookmarks")
          )
        ),
        fluidRow(
          bs4Dash::column(
            width = 8,
            offset = 2,
            br(), br(),
            tags$a(
              id = "start_report",
              class = "btn btn-default shiny-download-link",
              href = "",
              target = "_blank",
              download = NA,
              icon("download"),
              "Start the generation of the report"
            )
          )
        )
      )
    })

    output$ui_bookmarks <- renderUI({
      tagList(
        fluidRow(
          column(
            width = 6,
            bs4InfoBoxOutput("infobox_book_genes",
                             width = 6
            ),
            h5("Bookmarked genes"),
            DT::dataTableOutput("bookmarks_genes"),
            downloadButton("btn_export_genes", label = "", class = "biocdlbutton"),
            actionButton("btn_reset_genes",
                         label ="",
                         icon = icon("trash"),
                         style = .actionButtonStyle),

            bs4Dash::box(
              title = "Manually add genes to bookmarks",
              collapsible = TRUE,
              collapsed = TRUE,
              id = "box_bm_genes",
              width = 12,
              shinyAce::aceEditor(
                "editor_bookmarked_genes",
                theme = "solarized_light",
                height = "200px",
                readOnly = FALSE,
                wordWrap = TRUE,
                placeholder = paste(
                  c(
                    "Enter some gene identifiers or gene names, ",
                    "as they are provided in the `gene_id` and ",
                    "`gene_name` columns of the annotation object. ",
                    "For example:"
                  ),
                  collapse = "\n"
                )
              ),
              actionButton("load_bookmarked_genes",
                           label = "Upload genes",
                           icon = icon("upload"),
                           style = .actionButtonStyle)
            )
          ),
          column(
            width = 6,
            bs4InfoBoxOutput("infobox_book_genesets",
                             width = 6
            ),
            h5("Bookmarked genesets"),
            DT::dataTableOutput("bookmarks_genesets"),
            downloadButton("btn_export_genesets", label = "", class = "biocdlbutton"),
            actionButton("btn_reset_genesets",
                         label ="",
                         icon = icon("trash"),
                         style = .actionButtonStyle),

            bs4Dash::box(
              title = "Manually add genesets to bookmarks",
              collapsible = TRUE,
              collapsed = TRUE,
              id = "box_bm_genesets",
              width = 12,
              shinyAce::aceEditor(
                "editor_bookmarked_genesets",
                theme = "solarized_light",
                height = "200px",
                readOnly = FALSE,
                wordWrap = TRUE,
                placeholder = paste(
                  c(
                    "Enter some geneset identifiers, as they are ",
                    "provided in the `gs_id` column of the ",
                    "res_enrich object. You can also use the ",
                    "values in the `gs_description` field. ",
                    "For example:"
                  ),
                  collapse = "\n"
                )
              ),
              actionButton("load_bookmarked_genesets",
                           label = "Upload genesets",
                           icon = icon("upload"),
                           style = .actionButtonStyle)
            )
          )
        )
      )
    })

    output$infobox_book_genes <- renderbs4InfoBox({
      bs4InfoBox(
        title = "Bookmarked genes",
        value = length(reactive_values$genes),
        icon = icon("bookmark"),
        iconElevation = 3,
        color = "info",
        fill = TRUE,
        width = 12
      )
    })

    output$infobox_book_genesets <- renderbs4InfoBox({
      bs4InfoBox(
        title = "Bookmarked genesets",
        value = length(reactive_values$genes),
        icon = icon("bookmark"),
        iconElevation = 3,
        color = "success",
        fill = TRUE,
        width = 12
      )
    })

    output$bookmarks_genes <- DT::renderDataTable({
      DT::datatable(reactive_values$genesets, rownames = FALSE)
    })
    output$bookmarks_genesets <- DT::renderDataTable({
     DT::datatable(reactive_values$genesets, rownames = FALSE)
    })

    output$btn_export_genes <- downloadHandler(
      filename = function() {
        paste0("GeneTonicBookmarks_genes_", gsub(" ", "_", gsub("-", "", gsub(":", "-", as.character(Sys.time())))), ".txt")
      }, content = function(file) {
        writeLines(
          text = reactive_values$mygenes,
          con = file
        )
      }
    )

    output$btn_export_genesets <- downloadHandler(
      filename = function() {
        paste0("GeneTonicBookmarks_genesets_", gsub(" ", "_", gsub("-", "", gsub(":", "-", as.character(Sys.time())))), ".txt")
      }, content = function(file) {
        writeLines(
          text = reactive_values$mygenesets,
          con = file
        )
      }
    )

    output$start_report <- downloadHandler(
      filename = paste0(
        Sys.Date(),
        "_", round(stats::runif(1) * 100), # for not having all w the same name
        "_GeDiReport.html"
      ),
      content = function(file) {
        # temporarily switch to the temp dir, in case you do not have write permission to the current working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        withProgress(rmarkdown::render(
          input = system.file("extdata", "report_GeDi.Rmd", package = "GeDi"),
          output_file = file,
          quiet = TRUE
        ),
        message = "Generating the html report",
        detail = "This can take some time"
        )
      }
    )


    # Observers ------------------------------------------------------------------

    # Data Upload Panel ----------------------------------------------------------
    observeEvent(input$uploadgenesetfile, {
      filename <-
        unlist(strsplit(input$uploadgenesetfile$datapath, ".", fixed = TRUE))
      format <- tail(filename, n = 1)

      if (format == "RDS") {
        reactive_values$genesets <-
          readRDS(input$uploadgenesetfile$datapath)
      } else if (format == "txt") {
        reactive_values$genesets <- readGenesetsTxt()
      } else if (format == "xlsx") {
        reactive_values$genesets <-
          readxl::read_excel(input$uploadgenesetfile$datapath)
      } else {
        showNotification(
          "It seems like your input file has not the right format.
          Please check the Welcome panel for the right input format and
          upload your data again.",
          type = "error"
        )
      }

      reactive_values$gs_names <- reactive_values$genesets$Geneset
      tryCatch(
        expr = {
            reactive_values$genes <- getGenes(reactive_values$genesets)
        },
        error = function(cond) {
          showNotification(
            "It seems like your data does not have a column named 'Genes'. Please check your data and try to upload it again.",
            type = "error"
          )
          reactive_values$alt_names <- TRUE
          reactive_values$gs_names <- NULL
          reactive_values$genes <- NULL
        }
      )
    })

    observeEvent(input$alt_names_start, {
      reactive_values$gs_names <-
        reactive_values$genesets[, input$alt_name_genesets]

      reactive_values$genes <- getGenes(
        reactive_values$genesets,
        input$alt_name_genes
      )

      names(reactive_values$genesets)[names(reactive_values$genesets) == input$alt_name_genesets] <- "Genesets"
      names(reactive_values$genesets)[names(reactive_values$genesets) == input$alt_name_genes] <- "Genes"

      showNotification(
        "Successfully selected your columns and uploaded your data.
        We renamed the columns to 'Genesets' and 'Genes'.
        You can now proceed.",
        type = "message"
      )
    })

    observeEvent(input$btn_loaddemo, {
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Loading demo data now", value = 0)

      data(macrophage_topGO_example,
        package = "GeDi",
        envir = environment()
      )

      progress$inc(1 / 3, detail = "Extracting Genesets")

      reactive_values$genesets <- macrophage_topGO_example
      reactive_values$gs_names <- macrophage_topGO_example$Genesets

      progress$inc(1 / 3, detail = "Extracting Genes")
      reactive_values$genes <- getGenes(reactive_values$genesets)

      progress$inc(1 / 3, detail = "Successfully loaded demo data")
    })

    observeEvent(input$plot_brush, {
      df <- .buildHistogramData(
        genes = reactive_values$genes,
        gs_names = reactive_values$gs_names,
        start = input$bins_gs_hist[1],
        end = input$bins_gs_hist[2]
      )

      info_plot <- brushedPoints(df, input$plot_brush)
      output$table <- DT::renderDataTable(info_plot)
    })

    observeEvent(input$filter_genesets, {
        filtered_data <- .filterGenesets(
          input$select_filter_genesets,
          reactive_values$genesets
        )

      reactive_values$genesets <- filtered_data$Geneset
      reactive_values$gs_names <- filtered_data$gs_names
      reactive_values$genes <- filtered_data$Genes

      showNotification("Successfully filtered the selected Genesets.",
        type = "message"
      )

      # TODO: Does not work correctly currently, box still closes after one filtering round
      # updateBox(
      #   "optional_filtering_box",
      #   action = "update",
      #   options = list(
      #     id = "optional_filtering_box",
      #     width = 12,
      #     title = "Optional Filtering Step",
      #     status = "info",
      #     solidHeader = TRUE,
      #     h2("Filter your uploaded Genesets"),
      #     collapsible = TRUE,
      #     collapsed = FALSE
      #   )
      # )
    })

    observeEvent(input$download_ppi, {
      validate(need(!(input$species == ""),
        message = "Please specify the species of your data"
      ))
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
      id <- getId(reactive_values$species)
      progress$inc(1 / 12, detail = "Validate species ID")
      validate(
        need(
          !is.na(id),
          message = "We could not find your specified species.
                   Please check the spelling and try again."
        )
      )
      progress$inc(1 / 12, detail = "Get species specific STRINGdb")
      stringdb <- getStringDB(as.numeric(id))
      stringdb

      progress$inc(4 / 12, detail = "Get Annotation information")
      anno_df <- getAnnotation(stringdb)

      progress$inc(1 / 12, detail = "Download PPI")
      reactive_values$ppi <- getPPI(reactive_values$genes,
        string_db = stringdb,
        anno_df = anno_df
      )
      progress$inc(4 / 12, detail = "Successfully downloaded PPI")
    })

    # Scoring panel --------------------------------------------------------------
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
          "It seems like the file you've provided does not
             contain any genesets. Please check you input and retry.",
          type = "error"
        )
      }

      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Scoring your genesets", value = 0)

      if (input$scoringmethod == "Meet-Min") {
        scores <- getMeetMinMatrix(reactive_values$genes, progress)
      } else if (input$scoringmethod == "Kappa") {
        scores <- getKappaMatrix(reactive_values$genes, progress)
      } else if (input$scoringmethod == "PMM") {
        # TODO: Handle alpha as possible input
        if (is.null(reactive_values$ppi)) {
          showNotification(
            "It seems like you have not downloaded a PPI matrix. Please return to the Data Upload panel and download the respective PPI.",
            type = "error"
          )
          scores <- NULL
        } else {
          scores <- getpMMMatrix(reactive_values$genes,
            reactive_values$ppi,
            progress = progress
          )
        }
      } else if (input$scoringmethod == "Jaccard") {
        scores <- getJaccardMatrix(reactive_values$genes,
          progress = progress
        )
      }

      if (is.null(scores)) {
        showNotification(
          "It seems like something went wrong while scoring your data. Most likely this is due to the genesets not containg genes.",
          type = "error"
        )
      } else {
        rownames(scores) <- colnames(scores) <- reactive_values$gs_names
        reactive_values$scores <- scores
        updateBox("distance_calc_box", action = "toggle")
      }
    })

    observeEvent(input$create_heatmap,{
      InteractiveComplexHeatmap::InteractiveComplexHeatmapWidget(input,
                                                               output,
                                                               session,
                                                               scores_heatmap_react(),
                                                               output_id = "scores_heatmap",
                                                               width1 = 600,
                                                               height1 = 600,
                                                               width2 = 600,
                                                               height2 = 600,
                                                               click_action = .click_action,
                                                               brush_action = .brush_action,
                                                               close_button = FALSE,
                                                               output_ui = htmlOutput("info"))
    }
                 )

    .brush_action = function(df, output) {
      row = unique(unlist(df$row_index))
      column = unique(unlist(df$column_index))
      output[["info"]] = renderUI({
        if(!is.null(df)) {
          # HTML(kable_styling(kbl(as.matrix(reactive_values$scores)[row, column, drop = FALSE], digits = 2, format = "html"),
          #                    full_width = FALSE,
          #                    position = "left"))

          subset <- as.matrix(reactive_values$scores)[row, column, drop = FALSE]
          df <- as.data.frame(subset)
          DT::datatable(df,
                        options = list(scrollX = TRUE, scrollY = "400px"))
        }
      })
    }

    observeEvent(input$scores_graph_search, {
      current_node <- input$scores_graph_search
      if (current_node != "") {
        isolate({
          if (current_node %in% reactive_values$gs_names) {
            visNetworkProxy("scores_Network") %>% visSelectNodes(id = current_node)
          }
        })
      } else {
        showNotification(
          "It seems like the geneset identifier you searched for cannot be found
          in your input data. Please select a different geneset.",
          type = "message"
        )
      }
    })


    # Clustering panel -----------------------------------------------------------
    observeEvent(input$cluster_data, {
      if (is.null(reactive_values$scores)) {
        showNotification(
          "It seems like you did not score your Genesets yet. Please go back to
          the Scores panel and select a score of your choice.",
          type = "error"
        )
      }
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Finding initial seeds", value = 0)

      seeds <- seedFinding(
        reactive_values$scores,
        input$simThreshold,
        input$memThreshold
      )

      progress$inc(0.5, detail = "Found initial seeds. Now clustering the data.")

      cluster <- clustering(seeds, input$clustThreshold)

      progress$inc(0.5, detail = "Finished clustering the data.")

      if (is.null(seeds)) {
        showNotification(
          "It seems like something went wrong while clustering your data. Most likely this is due to the genesets not containg genes.",
          type = "error"
        )
      } else {
        reactive_values$seeds <- seeds
        reactive_values$cluster <- cluster
        bs4Dash::updateBox("clustering_param_box", action = "toggle")
      }
    })


    # Tour Observers -------------------------------------------------------------
    observeEvent(input$tour_welcome, {
      tour <- utils::read.delim(
        system.file("extdata",
          "intro_welcome.txt",
          package = "GeDi"
        ),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_data_upload, {
      tour <- utils::read.delim(
        system.file("extdata",
          "intro_data_upload.txt",
          package = "GeDi"
        ),
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
          package = "GeDi"
        ),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour, tooltipPosition = "bottom"))
    })

    observeEvent(input$tour_graph, {
      tour <- utils::read.delim(
        system.file("extdata",
                    "intro_graph.txt",
                    package = "GeDi"
        ),
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })


    #Buttons Navbar Observers -------------------------------------------------
    observeEvent(input$btn_first_help, {
      showModal(
        modalDialog(
          title = "First Help Info", size = "l", fade = TRUE,
          footer = NULL, easyClose = TRUE,
          tagList(
            includeMarkdown(system.file("extdata", "GeDi101.md", package = "GeDi")),
          )
        )
      )
    })

    observeEvent(input$btn_info_session, {
      showModal(
        modalDialog(
          title = "Session information", size = "l", fade = TRUE,
          footer = NULL, easyClose = TRUE,
          tagList(
            tags$code("> sessionInfo()"),
            renderPrint({
              utils::sessionInfo()
            })
          )
        )
      )
    })

    observeEvent(input$btn_info_gedi, {
      showModal(
        modalDialog(
          title = "About GeDi", size = "l", fade = TRUE,
          footer = NULL, easyClose = TRUE,
          tagList(
            includeMarkdown(system.file("extdata", "about.md", package = "GeDi")),
            renderPrint({
              utils::citation("GeneTonic")
            })
          )
        )
      )
    })

  }
  shinyApp(ui = gedi_ui, server = gedi_server)
}
