#' Create a minimal chevreulShiny app using SingleCellExperiment input
#'
#' @param single_cell_sce a singlecell object
#' @param appTitle a title for the app
#' @param organism_type human or mouse
#' @param futureMb the megabytes available for the future package
#' @param db_name a database of bigwig files
#'
#' @return a minimal chevreulShiny app
#' @export
#'
#' @examples
#' 
#' if (interactive() ) {
#'   data("tiny_sce")
#'   minimalChevreulApp(tiny_sce)
#' }
#' 
minimalChevreulApp <- function(single_cell_sce = NULL,
                          appTitle = NULL,
                          organism_type = "human",
                          futureMb = 13000,
                          db_name = "single-cell-projects.db",
                          bigwig_db = bw_files) {
    
    db_path <- file.path(user_cache_dir(appname="chevreul"), db_name)

    get_adjusted_cores <- function() {
        available_cores <- availableCores()[1]
        if (available_cores > 2) {
            return(available_cores - 2)
        } else {
            return(1)
        }
    }
    
    if(.Platform$OS.type == "windows"){
        plan(strategy = "multisession", workers = get_adjusted_cores())  # Limit to 6 cores
    } else {
        plan(strategy = "multicore", workers = get_adjusted_cores())  # Limit to 6 cores
    }
    
    # sets the maximum allowed total size (in bytes) of global variables
    future_size <- futureMb * 1024^2
    options(future.globals.maxSize = future_size)
    
    # controls the maximum file size allowed for uploads
    options(shiny.maxRequestSize = 40 * 1024^2)

    options(DT.options = list(
        pageLength = 2000, paging = FALSE,
        info = TRUE, searching = TRUE, autoWidth = FALSE, ordering = TRUE,
        scrollX = TRUE, language = list(search = "Filter:")
    ))
    header <- dashboardHeader(title = appTitle)
    sidebar <- dashboardSidebar(
        textOutput("appTitle"),
        uiOutput("featureType"),
        sidebarMenu(
            menuItem("Overview Plots",
                tabName = "comparePlots", icon = icon("chart-bar"),
                selected = TRUE
            ), menuItem("Violin",
                tabName = "violinPlots", icon = icon("sort")
                ), menuItem("Differential Expression",
                tabName = "diffex", icon = icon("magnet")
            ), menuItem("Find Markers",
                tabName = "findMarkers", icon = icon("bullhorn")
            ), menuItem("Heatmap",
                       tabName = "heatmapPlots", icon = icon("square")
            ), menuItem("Transcripts",
                        tabName = "allTranscripts", icon = icon("sliders-h")
            ), menuItem("Subset",
                tabName = "subsetObject", icon = icon("filter")
            # ),  menuItem("Regression",
            #     tabName = "regressFeatures", icon = icon("eraser")
            ),menuItem("Formatting",
                       tabName = "reformatMetadata", icon = icon("columns")
            ), menuItem("Coverage Plots",
                        tabName = "coveragePlots", icon = icon("mountain")
            ), menuItem("Technical Information",
                tabName = "techInfo", icon = icon("cogs")
            )
        ),
        actionButton("changeEmbedAction",
            label = "Change Embedding Parameters"
        ),
        changeEmbedParamsui("changeembed"),
        width = 250
    )
    body <- dashboardBody(
        use_waiter(),
        tabItems(
            tabItem(
                tabName = "comparePlots",
                h2("Overview Plots") |>
                    default_helper(type = "markdown", 
                                   content = "comparePlots"),
                plotDimRedui("plotdimred1"),
                plotDimRedui("plotdimred2"),
                plotReadCountui("plotreadcount1"),
                plotReadCountui("plotreadcount2"),
                chevreulBox(
                    title = "Selected Cells",
                    tableSelectedui("tableselected"),
                    width = 6
                ),
                plotClustree_UI("clustreePlot")
            ),
            tabItem(
                tabName = "violinPlots",
                fluidRow(
                    plotViolinui("violinPlot")
                )
            ),
            tabItem(
                tabName = "heatmapPlots",
                fluidRow(
                    plotHeatmapui("heatMap")
                )
            ),
            tabItem(
                tabName = "coveragePlots",
                fluidRow(
                    plotCoverage_UI("coverageplots")
                )
            ),
            tabItem(
                tabName = "reformatMetadata",
                fluidRow(
                    reformatMetadataDRui("reformatMetadataDR")
                )
            ),
            tabItem(
                tabName = "subsetObject",
                h2("Subset Object Input") |>
                    default_helper(type = "markdown", 
                                   content = "subsetObject"),
                plotDimRedui("subset"),
                chevreulBox(
                    title = "Subset Settings",
                    actionButton("resetObject",
                                 "reset object to initial state"),
                    actionButton("subsetAction",
                                 "subset object by selected cells"),
                    actionButton("subsetCsv", "subset object by uploaded csv"),
                    fileInput("uploadCsv",
                        "Upload .csv file with cells to include",
                        accept = c(".csv")
                    ),
                    useShinyjs(),
                    # runcodeUI(code = "alert('Hello!')"),
                    textOutput("subsetMessages"),
                    width = 6
                ),
                chevreulBox(
                    title = "Selected Cells", tableSelectedui("subset"),
                    width = 6
                )
            ), tabItem(
                tabName = "findMarkers",
                h2("Find Markers"),
                chevreulMarkersui("findmarkers"),
                plotDimRedui("markerScatter")
            ), tabItem(
                tabName = "allTranscripts",
                h2("All Transcripts"),
                plotDimRedui("alltranscripts2"),
                allTranscriptsui("alltranscripts1")
            ),
            tabItem(
                tabName = "diffex",
                h2("Differential Expression") |>
                    default_helper(type = "markdown", content = "diffex"),
                plotDimRedui("diffex"),
                diffexui("diffex")
            ),
            tabItem(
                tabName = "regressFeatures",
                fluidRow(
                    chevreulBox(
                        title = "Regress Features",
                        actionButton(
                            "regressAction",
                            "Regress Objects By Genes"
                        ),
                        width = 12
                    ) |>
                        default_helper(type = "markdown",
                                       content = "regressFeatures")
                )),
            tabItem(
                tabName = "techInfo",
                h2("Technical Information"),
                techInfoui("techInfo")
            )
        )
    )

    ui <- function(request) {
        ui <- dashboardPage(
            header = header, sidebar = sidebar,
            body = body
        )
    }
    server <- function(input, output, session) {
        w <- Waiter$new()

        observe_helpers(
            help_dir = system.file("helpers", package = "chevreulShiny"))
        options(warn = -1)

        object <- reactiveVal(NULL)
        observe({
            req(!is.null(single_cell_sce))
            object(single_cell_sce)
        })

        organism_type <- reactive({
            "human"
        })

        plot_types <- reactive({
            req(!is.null(object()))
            list_plot_types(object())
        })

        featureType <- reactive({
            "gene"
            # input$feature_type
        })


        observe({
            reformatted_sce <- reformatMetadataDR( 
                                             "reformatMetadataDR", 
                                             object, featureType)
            object(reformatted_sce())
        })

        reductions <- reactive({
            req(!is.null(object()))
            reducedDimNames(object())
        })

        observe({
            req(!is.null(object()))

            plotDimRed("plotdimred1", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
            plotDimRed("plotdimred2", object, plot_types, 
                       featureType,
                organism_type = organism_type, reductions
            )
            plotDimRed("diffex", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
            plotDimRed("subset", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
            plotDimRed("markerScatter", object, plot_types, 
                       featureType,
                organism_type = organism_type, reductions
            )
            plotReadCount("plotreadcount1", object, plot_types)
            plotReadCount("plotreadcount2", object, plot_types)
        })
        plotViolin("violinPlot", object, featureType,
            organism_type
        )
        plotHeatmap("heatMap", object, featureType,
            organism_type
        )
        plotCoverage("coverageplots", object, plot_types, proj_dir, 
            organism_type, bigwig_db = bigwig_db
        )
        plotClustree("clustreePlot", object)
        tableSelected("tableselected", object)
        diffex_selected_cells <- tableSelected("diffex",
            object
        )
        subset_selected_cells <- tableSelected("subset",
            object
        )

        observeEvent(input$changeEmbedAction, {
            showModal(modalDialog(
                title = "Recalculating Embedding",
                "This process may take a minute or two!"
            ))
            object <- changeEmbedParams( "changeembed",
                object
            )
            removeModal()
        })

        chevreulMarkers("findmarkers", object, plot_types, 
                   featureType)

        # plot all transcripts
        observe({
            req(featureType())
            req(!is.null(object()))
            allTranscripts( "alltranscripts1", object, featureType,
                organism_type
            )

            plotDimRed( "alltranscripts2", object, plot_types, 
                       featureType,
                organism_type = organism_type, reductions
            )
        })

        diffex("diffex", object, featureType,
            diffex_selected_cells
        )
        
        # subset by lasso selection
        observeEvent(input$subsetAction, {
            req(input$subsetAction)
            req(subset_selected_cells())
            withCallingHandlers(
                {
                    html("subsetMessages", "")
                    message("Beginning")
                    subset_sce <- 
                        object()[, colnames(object()) %in% 
                                     subset_selected_cells()]
                    object(subset_sce)
                    if (length(unique(object()$batch)) > 1) {
                        message("reintegrating gene expression")
                        reintegrated_sce <- reintegrate_sce(object(),
                            resolution = seq(0.2, 1, by = 0.2),
                            organism = metadata(object())$experiment$organism
                        )
                        object(reintegrated_sce)
                    } else {
                        processed_sce <- 
                            sce_process(object(), 
                                            resolution = seq(0.2, 1, by = 0.2))
                        object(processed_sce)
                    }
                    message("Complete!")
                },
                message = function(m) {
                    html(id = "subsetMessages", html = paste0(
                        "Subsetting Object: ",
                        m$message
                    ), add = FALSE)
                }
            )
        })
        
        # subset by .csv
        observeEvent(input$subsetCsv, {
            req(input$subsetCsv)
            req(input$uploadCsv)
            withCallingHandlers(
                {
                    html("subsetMessages", "")
                    message("Beginning")
                    subset_sce <- subset_by_colData(
                        input$uploadCsv$datapath,
                        object()
                    )
                    object(subset_sce)

                    if (length(unique(object()[["batch"]])) > 1) {
                        message("reintegrating gene expression")
                        reintegrated_sce <- reintegrate_sce(object(),
                            resolution = seq(0.2, 1, by = 0.2),
                            organism = metadata(object())$experiment$organism
                        )
                        object(reintegrated_sce)
                    } else {
                        processed_sce <- 
                            sce_process(object(), 
                                           resolution = seq(0.2, 1, by = 0.2))
                        object(processed_sce)
                    }
                    message("Complete!")
                },
                message = function(m) {
                    html(id = "subsetMessages", html = paste0(
                        "Subsetting Object: ",
                        m$message
                    ), add = FALSE)
                }
            )
        })
        
        # reset to original object
        observeEvent(input$resetObject, {
            object(single_cell_sce)
        })

        observeEvent(input$regressAction, {
            req(!is.null(object()))
            showModal(modalDialog(
                title = "Regressing out cell cycle effects",
                "This process may take a minute or two!"
            ))
            regressed_sce <- regress_cell_cycle(object())
            object(regressed_sce)
            removeModal()
        })

        techInfo("techInfo", object)
    }
    shinyApp(ui, server, enableBookmarking = "server")
}
