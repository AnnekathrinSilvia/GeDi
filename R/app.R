#' GeDi main function
#'
#' @param genesets a dataframe of genesets containing the columns: genesets (the names / ids of the sets) and genes (the genes included in the genesets)
#' @param ppi a Protein-Protein interaction matrix
#' @param alpha a scaling factor in between 0 and 1
#'
#' @return A Shiny app object is returned
#' @export
#' @import DT
#' @import visNetwork
#' @import shiny
#' @import shinyBS
#' @import fontawesome
#' @importFrom rintrojs introjs
#' @importFrom utils read.delim
#' @importFrom bs4Dash bs4DashPage bs4DashNavbar box bs4DashBrand bs4DashBody bs4DashSidebar bs4SidebarMenu bs4SidebarMenuItem bs4TabItem bs4TabItems
#'
#' @examples
#' \dontrun{
#' GeDi()
#' }
GeDi <- function(genesets = NULL,
                 ppi = NULL,
                 alpha = 1) {
  oopt <- options(spinner.type = 6)
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
            label = "First Help"
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
            inputId = "btn_info_session",
            icon = icon("circle-info"),
            label = "About this session"
          ),
          actionButton(
            inputId = "btn_info_genetonic",
            icon = icon("heart"),
            label = "About GeDi"
          )
        )
      )
    ),
    # ),
    # title = bs4DashBrand(
    #   title = HTML("<small>GeDi</small>"),
    #   #href = "https://bioconductor.org/packages/GeneTonic",
    #   # color = "info",
    #   #image = "GeneTonic/GeneTonic.png"),
    #   skin = "dark",
    #   status = "gray-dark",
    #   border = FALSE,
    #   controlbarIcon = icon("cogs"),
    #   fixed = TRUE
    # ),

    # sidebar definition ------------------------------------------------------
    sidebar = bs4DashSidebar(
      id = "sidebar",
      title = HTML("<small>GeDi</small>"),
      #src = "GeneTonic/GeneTonic.png",
      skin = "dark",
      status = "primary",
      brandColor = NULL,
      #url = "https://bioconductor.org/packages/GeneTonic",
      collapsed = TRUE,
      elevation = 1,
      opacity = 0.8,
      bs4SidebarMenu(
        bs4SidebarMenuItem(
          "Welcome!",
          tabName = "tab_welcome",
          icon = icon("house")),
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
        bs4TabItem(tabName = "tab_welcome",
                   tagList(
                     fluidRow(column(width = 11),
                              column(
                                width = 1,
                                actionButton(
                                  "tour_welcome",
                                  label = "",
                                  icon = icon("circle-question"),
                                  style = "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
                                )
                              )),
                     uiOutput("ui_panel_welcome")
                   )),

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
                  style = "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
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
                           style = "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
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
                           style = "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
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
                   ))

      )
    )
    # ),
    # # controlbar definition ------------------------------------------------
    # controlbar = bs4Dash::bs4DashControlbar(
    #   collapsed = TRUE,
    #   uiOutput("ui_controlbar")
    # ),
    #
    # # footer definition -------------------------------------------------------
    # footer = bs4DashFooter(
    #   left = GeneTonic_footer,
    #   right = NULL
    # )
  )


  # server deifnition ---------------------------------------------------------
  gedi_server <- function(input, output, session) {
    # initializing reactives --------------------------------------------------
    reactive_values <- reactiveValues()

    # input data
    reactive_values$genesets <- NULL
    #only the names of the genesets -> TODO: evaluate later if needed
    reactive_values$gs_names <- NULL
    #splitted genes
    reactive_values$genes <- NULL
    reactive_values$species <- NULL
    reactive_values$ppi <- NULL
    reactive_values$scores <- NULL


    # panel Welcome ----------------------------------------------------------
    output$ui_welcome <- renderUI({
      tagList(fluidRow(column(width = 12)))
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
              uiOutput("upload_genesets"),
              br(),
              "... or you can also ",
              # actionButton("btn_loaddemo", "Load the demo airway data",
              #              icon = icon("play-circle"),
              #              class = "btn btn-info"
              # ), br(), p()),
            ),
            column(
              width = 6,
              box(
                id = "Genests_preview",
                width = NULL,
                title = "Genesets preview",
                status = "primary",
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

    output$upload_genesets <- renderUI({
      return(
        fileInput(
          inputId = "uploadgenesetfile",
          label = "Upload a geneset file",
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            "text/tab-separated-values",
            "text/plain",
            ".csv",
            ".tsv"
          ),
          multiple = FALSE
        )
      )
    })

    readGenesets <- reactive({
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

    output$dt_genesets <- DT::renderDataTable({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      DT::datatable(reactive_values$genesets,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })


    output$ui_panel_specify_species <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      box(
        width = 12,
        title = "Step 2",
        status = "warning",
        solidHeader = TRUE,
        tagList(
          h2("Select the species of your data"),
          fluidRow(column(
            width = 6,
            uiOutput("ui_specify_species")
          ),
          column(
            width = 6,
            uiOutput("ui_specify_opt_params_ppi_download")
          ))
        )
      )
    })

    output$ui_specify_species <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      fluidRow(column(
        width = 12,
        textInput(
          "species",
          label = "Please specify the species of your data
                          (e.g. Homo Sapiens, Mus musculus, etc.)",
          width = '400px'
        )
      ))
    })

    output$ui_specify_opt_params_ppi_download <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      box(
        id = "Optional_parameters",
        width = 12,
        title = "Optional Parameters",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        tagList(fluidRow(column(
          width = 6,
          uiOutput("ui_optional_param")
        )))
      )
    })

    output$ui_optional_param <- renderUI({
      fluidRow(column(
        width = 12,
        textInput("stringVersion",
                  label = "Specify the StringDb version to use",
                  value = "11.5"),
        textInput(
          "scoreThresholdString",
          label = "Specify the score threshold",
          value = "0.00"
        ),
        textInput(
          "inputDirectoryString",
          label = "Specify the input directory",
          value = ""
        )
      ))
    })

    output$ui_panel_download_ppi <- renderUI({
      if (input$species == "" ||
          is.na(input$species) || is.null(input$species)) {
        return(NULL)
      }
      reactive_values$species <- input$species
      box(
        width = 12,
        title = "Step 3",
        status = "success",
        solidHeader = TRUE,
        tagList(
          h2("Download the PPI matrix from STRING"),
          fluidRow(column(
            width = 6,
            actionButton(
              "download_ppi",
              label = "Download PPI matrix",
              icon = icon("download"),
              style = "color: #FFFFFF; background-color: #0092AC; border-color: #0092AC"
            )
          ),
          column(
            width = 6,
            box(
              id = "PPI_preview",
              width = NULL,
              title = "PPI preview",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              collapsed = TRUE,
              fluidRow(column(
                width = 12,
                offset = 0.5,
                DT::dataTableOutput("dt_ppi")
              ))
            )
          ))
        )
      )
    })

    output$dt_ppi <- DT::renderDataTable({
      if (is.null(reactive_values$ppi)) {
        return(NULL)
      }
      DT::datatable(reactive_values$ppi,
                    options = list(scrollX = TRUE, scrollY = "400px"))
    })

    # panel Scores -----------------------------------------------------

    output$ui_panel_scores <- renderUI({
      tagList(fluidRow(
        column(
          width = 11,
          selectInput(
            "scoringmethod",
            label = "Select the scoring method for your data: ",
            choices = c("", "Meet-Min", "Kappa", "PMM"),
            multiple = FALSE
          ),
          uiOutput("ui_score_data")
        ),
        column(width = 1)
      ),
      fluidRow(column(
        width = 11,
        visNetworkOutput("ui_plot_scores")
      )))
    })

    # TODO: Workaround with empty string, check to fix that later in another way
    output$ui_score_data <- renderUI({
      if (is.null(input$scoringmethod) || input$scoringmethod == "") {
        return(NULL)
      }
      fluidRow(
        column(
          width = 5,
          strong("Now you can score your data"),
          br(),
          "Attention: If you have many genesets to score,
           this operation may take some time",
          br(),
          actionButton("score_data",
                       label = "Score the Genesets",
                       style = "color: #FFFFFF; background-color: #0092AC; border-color: #0092AC")
        )
      )
    })

    output$ui_plot_scores <- renderVisNetwork({
      if (is.null(reactive_values$scores)) {
        return(NULL)
      }
      #pheatmap::pheatmap(reactive_values$scores)
      adj <- getAdjacencyMatrix(reactive_values$scores, 0.3)
      g <- buildGraph(adj)

      visNetwork::visIgraph(g) %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = TRUE
        ) %>%
        visExport(name = "backbone_network",
                  type = "png",
                  label = "Save backbone graph")
    })


    # panel Graph ------------------------------------------------------

    output$ui_panel_graph <- renderUI({
      tagList(fluidRow(column(width = 11),
                       column(width = 1)),
              fluidRow(column(width = 8,
                              fluidRow(
                                column(width = 1)
                              )),
                       column(width = 4)))
    })

    # controlbar --------------------------------------------------------------

    # output$ui_controlbar <- renderUI({
    #   validate(
    #     need(!is.null(reactive_values$gtl) ,
    #          message = "\n\n\nPlease provide a GeneTonicList or its components to display the content of the control bar.")
    #   )
    #
    #   # message(nrow(reactive_values$res_enrich))
    #
    #   tagList(
    #     numericInput(
    #       inputId = "de_fdr",
    #       label = "False Discovery Rate (FDR) for DE",
    #       value = 0.05, min = 0.0001, max = 1, step = 0.01
    #     ),
    #     numericInput(
    #       inputId = "n_genesets",
    #       label = "Number of genesets",
    #       value = 15, min = 1, max = nrow(reactive_values$res_enrich)
    #     ),
    #     selectInput("exp_condition",
    #                 label = "Group/color by: ",
    #                 choices = c(NULL, names(colData(reactive_values$dds))), selected = NULL, multiple = TRUE
    #     ),
    #     colourInput("col", "Select colour for volcano plot", "#1a81c2",
    #                 returnName = TRUE,
    #                 allowTransparent = TRUE
    #     ),
    #     checkboxInput("labels", label = "Display all labels", value = FALSE)
    #   )
    # })

    # outputOptions(output, "ui_controlbar", suspendWhenHidden = FALSE)

    # Observers ------------------------------------------------------------------------------------
    observeEvent(input$uploadgenesetfile, {
      reactive_values$genesets <- readGenesets()
      reactive_values$gs_names <- reactive_values$genesets$Geneset
      tryCatch(
        expr = {
          reactive_values$genes <- getGenes(reactive_values$genesets)
        },
        error = function(cond) {
          message("It seems like your data does not have a column named 'Genes'.")
          message("Here's the original error message:")
          message(cond)
          message("Please check your data and try to upload it again.")
          reactive_values$genesets <- NULL
          reactive_values$gs_names <- NULL
          reactive_values$genes <- NULL
        }
      )
    })

    observeEvent(input$download_ppi, {
      validate(need(!(input$species == ""),
                    message = "Please specify the species of your data"))
      reactive_values$species <- input$species
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Downloading PPI from STRINGdb", value = 0)
      progress$inc(1/12, detail = paste("Get species ID"))
      id <- getId(reactive_values$species)
      progress$inc(1/12, detail = paste("Validate species ID"))
      validate(
        need(
          !is.na(id),
          message = "We could not find your specified species.
                   Please check the spelling and try again."
        )
      )
      progress$inc(1/12, detail = paste("Get species specific STRINGdb"))
      stringdb <-
        getStringDB(as.numeric(id))
      stringdb

      progress$inc(4/12, detail = paste("Get Annotation information"))
      anno_df <- getAnnotation(stringdb)

      progress$inc(1/12, detail = paste("Download PPI"))
      reactive_values$ppi <-
        getPPI(reactive_values$genes,
               string_db = stringdb,
               anno_df = anno_df)

      progress$inc(4/12, detail = paste("Done"))
    })

    observeEvent(input$score_data, {
      validate(need(
        input$scoringmethod != "" &
          !(is.null(input$scoringmethod)),
        "Please select a scoring method"
      ))
      validate(
        need((length(reactive_values$genes) > 0  &&
                !(is.null(
                  reactive_values$genes
                ))) ,
             message = "\n\n\nIt seems like the file you've provided does not
             contain any genesets. Please check you input and retry."
        )
      )
      if (input$scoringmethod == "Meet-Min") {
        scores <-
          getMeetMinMatrix(reactive_values$genes)
      } else if (input$scoringmethod == "Kappa") {
        scores <-
          getKappaMatrix(reactive_values$genes)
      } else if (input$scoringmethod == "PMM") {
        # TODO: Handle alpha as possible input
        scores <-
          getpMMMatrix(reactive_values$genes,
                       reactive_values$ppi)
      }

      validate(
        need((scores != -1) ,
             message = "\n\n\nIt seems like the file you've provided does not
             contain any genesets. Please check you input and retry."
        )
      )

      rownames(scores) <-
        colnames(scores) <- reactive_values$gs_names
      reactive_values$scores <- scores
    })


    # Tour Observers -------------------------------------------------------------------------
    observeEvent(input$tour_welcome, {
      tour <- utils::read.delim(
        "tours/intro_welcome.txt",
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_data_upload, {
      tour <- utils::read.delim(
        "tours/intro_data_upload.txt",
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_scoring, {
      tour <- utils::read.delim(
        "tours/intro_scoring.txt",
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

    observeEvent(input$tour_graph, {
      tour <- utils::read.delim(
        "tours/intro_graph.txt",
        sep = ";",
        stringsAsFactors = FALSE,
        row.names = NULL,
        quote = ""
      )
      introjs(session, options = list(steps = tour))
    })

  }
  shinyApp(ui = gedi_ui, server = gedi_server)
}
