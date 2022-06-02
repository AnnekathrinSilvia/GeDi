#' GeDi main function
#'
#' @param genesets a dataframe of genesets containing the columns: genesets (the names / ids of the sets) and genes (the genes included in the genesets)
#' @param ppi a Protein-Protein interaction matrix
#' @param alpha a scaling factor in between 0 and 1
#'
#' @return
#' @export
#' @import DT
#' @import igraph
#' @import visNetwork
#' @import shiny
#' @import shinyBS
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
            icon = icon("question-circle"),
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
            icon = icon("info-circle"),
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
        bs4SidebarMenuItem("Welcome!",
                           tabName = "tab_welcome",
                           icon = icon("home")),
        bs4SidebarMenuItem(
          "Data Upload",
          tabName = "tab_data_upload",
          icon = icon("share-alt-square")
        ),
        bs4SidebarMenuItem(
          "Distance Scores",
          tabName = "tab_scores",
          icon = icon("project-diagram")
        ),
        bs4SidebarMenuItem(
          "Graph",
          tabName = "tab_graph",
          icon = icon("share-alt-square")
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
            "
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
                   uiOutput("ui_panel_welcome")),

        # ui panel data upload -----------------------------------------------
        bs4TabItem(
          tabName = "tab_data_upload",
          uiOutput("ui_panel_data_upload"),
          uiOutput("ui_download_ppi")
        ),

        # ui panel scores -------------------------------------------------
        bs4TabItem(tabName = "tab_scores",
                   uiOutput("ui_panel_scores")),

        # ui panel graph ---------------------------------------------------
        bs4TabItem(tabName = "tab_graph",
                   uiOutput("ui_panel_graph"))

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
    # original input data
    reactive_values$genesets <- NULL
    #only the names of the genesets -> TODO: evaluate later if needed
    reactive_values$gs_names <- NULL
    #splitted genes
    reactive_values$genes <- NULL
    reactive_values$genes_ind <- NULL
    reactive_values$ppi <- NULL
    reactive_values$scores <- NULL


    # panel Welcome ----------------------------------------------------------
    output$ui_welcome <- renderUI({
      tagList(fluidRow(column(width = 12)))
    })


    # panel Data Upload ------------------------------------------------------

    output$ui_panel_data_upload <- renderUI({
      tagList(
        fluidRow(column(
          width = 8,
          bsCollapse(
            id = "help_datasetup",
            open = NULL,
            bsCollapsePanel("Help")
          )
        )),
        box(
          width = 12,
          title = "Step 1",
          status = "danger",
          solidHeader = TRUE,
          h2("Upload your Genesets input data"),

          fluidRow(
            column(
              width = 4,
              uiOutput("upload_genesets"),
              br(),
              "... or you can also ",
              # actionButton("btn_loaddemo", "Load the demo airway data",
              #              icon = icon("play-circle"),
              #              class = "btn btn-info"
              # ), br(), p()),
              column(
                width = 4,
                br(),
                actionButton(
                  "help_format",
                  label = "",
                  icon = icon("question-circle"),
                  style = "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
                ),
                bsTooltip(
                  "help_format",
                  "How to provide your input data to ideal",
                  "bottom",
                  options = list(container = "body")
                )
              )
            ),

            fluidRow(column(
              width = 12,
              box(
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
            ))
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
      # TODO: Build own sepguesser so ideal dependency can be removed
      guessed_sep <-
        ideal::sepguesser(input$uploadgenesetfile$datapath)
      genesets <-
        utils::read.delim(
          input$uploadgenesetfile$datapath,
          header = TRUE,
          as.is = TRUE,
          sep = guessed_sep,
          quote = "",
          row.names = 1,
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


    output$ui_download_ppi <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      box(
        width = 12,
        title = "Step 2",
        status = "warning",
        solidHeader = TRUE,
        tagList(
          h2(
            "Select the species of you data and download the corresponding
             Protein-Protein Interaction (PPI) Matrix"
          ),
          fluidRow(
            column(width = 6,
                   uiOutput("ui_species")),
            column(
              width = 6,
              strong("Download the PPI matrix from STRING:"),
              br(),
              "Attention: This operation can take some time",
              br(),
              actionButton(
                "download_ppi",
                label = "Download PPI matrix",
                icon = icon("download"),
                style = "color: #FFFFFF; background-color: #0092AC; border-color: #0092AC"

              )
            )

          ),
          fluidRow(column(
            width = 6,
            uiOutput("ui_opt_param_ppi")
          ))
        )
      )
    })

    output$ui_species <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      poss_species <-
        read.delim("data/species.txt", sep = "\n", header = F)
      selectInput(
        "species",
        label = "Select the species of your data: ",
        choices = poss_species,
        selected = poss_species[1],
        multiple = FALSE
      )
    })

    output$ui_opt_param_ppi <- renderUI({
      if (is.null(reactive_values$genesets)) {
        return(NULL)
      }
      box(
        width = 12,
        title = "Optional Parameters",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        tagList(fluidRow(column(
          width = 6,
          uiOutput("ui_cachePath")
        )))
      )
    })

    output$ui_cachePath <- renderUI({
      textInput("cachePath",
                label = "Change the path for the PPI Cache",
                value = "~GSD",
                width = '400px')

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

    # TODO: Workaround with empty string check to fiy that later in another way
    output$ui_score_data <- renderUI({
      if (is.null(input$scoringmethod) || input$scoringmethod == "") {
        return(NULL)
      }
      fluidRow(
        column(
          width = 5,
          strong("Now you can score your data"),
          br(),
          "Attention: If you have many genesets to score, this operation may take some time",
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
        visExport(
          name = "backbone_network",
          type = "png",
          label = "Save backbone graph"
        )
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


    observeEvent(input$uploadgenesetfile, {
      reactive_values$genesets <- readGenesets()
      reactive_values$gs_names <- rownames(reactive_values$genesets)
      reactive_values$genes <- getGenes(reactive_values$genesets)
      #reactive_values$genes_ind <- getGenesIndexed(reactive_values$genes)
    })

    observeEvent(input$download_ppi, {
      if (is.null(input$cachePath)) {
        reactive_values$ppi <- downloadPPI(input$species)
      } else{
        reactive_values$ppi <-
          downloadPPI(input$species, cachepath = input$cachePath)
      }
    })

    observeEvent(input$score_data, {
      validate(need(
        input$scoringmethod != "" &
          !(is.null(input$scoringmethod)),
        "Please select a scoring method"
      ))
      if (input$scoringmethod == "Meet-Min") {
        reactive_values$scores <- getMeetMinMatrix(reactive_values$genes, reactive_values$gs_names)
      } else if (input$scoringmethod == "Kappa") {
        reactive_values$scores <- getKappaMatrix(reactive_values$genes, reactive_values$gs_names)
      } else if (input$scoringmethod == "PMM") {
        reactive_values$scores <-
          getpMMMatrix(reactive_values$genes, reactive_values$ppi, reactive_values$gs_names)
      }

      validate(
        need((reactive_values$scores != -1) ,
             message = "\n\n\nIt seems like the file you've provided does not contain any genesets. Please check you input and retry.")
       )
    })

  }
  shinyApp(ui = gedi_ui, server = gedi_server)
}
