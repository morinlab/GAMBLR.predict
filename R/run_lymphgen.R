#' @title Classify DLBCL samples with LymphGen
#'
#' @description Run LymphGen (Wright et al., Cancer Cell 2020) on the inputs
#' built by [lymphgen_input], with a local copy of the LymphGen R
#' distribution and its training data. LymphGen is not part of GAMBLR: obtain
#' it separately and point `lymphgen_dir` at it.
#'
#' @details `lymphgen_dir` must hold LymphGen 2's command-line distribution:
#' `Run_Lymphgen.R`, `Source LymphGen.R` and its training and reference files
#' (Fullmat.2020.txt, Full.Uber.2020.txt, Featord.2020.txt, Index.Flat.txt,
#' Study.samp.txt, originalpred.txt / originalpredv2.txt and Result.txt).
#' Give it to each call, or set it for the session with
#' `options(GAMBLR.predict.lymphgen_dir = ...)` or the environment variable
#' `LYMPHGEN_DIR`. `Run_Lymphgen.R` needs the optparse package in the R that
#' runs it (`rscript`).
#'
#' The input files are written to a temporary directory and classified by
#' `Run_Lymphgen.R`, so the calls are LymphGen's. Without copy number
#' (`copy_flat` NULL), LymphGen runs its no-copy-number model. `Run_Lymphgen.R`
#' doesn't write out the warnings it collects (e.g. genes in the flat file
#' but not in the gene list); its console output is in
#' `attr(result, "lymphgen_run")$log`.
#'
#' @param lymphgen_input The output of [lymphgen_input] (its mutation_flat,
#' gene_list, sample_annotation and, with copy number, copy_flat and arm), or
#' a list of them (unnamed ones are named input1, input2, ...).
#' @param lymphgen_dir Directory of the LymphGen distribution. Default: the
#' `GAMBLR.predict.lymphgen_dir` option, else `LYMPHGEN_DIR`.
#' @param model LymphGen model: "original" (LymphGen 2 as released) or
#' "gambl_v2", the GAMBL-derived model (no copy number, no A53, composite
#' calls at 0.40; `Run_Lymphgen.R --model gambl_v2`), which needs the GAMBL
#' LymphGen repository (with `models/gambl_v2/`) as `lymphgen_dir`. Default
#' NULL: the model the inputs were built for (attribute "lymphgen_model" of
#' [lymphgen_input]'s output), else "original". For gambl_v2 the result also
#' has LymphGen.Default.Prediction, LymphGen's own call. "gambl_v3" extends
#' gambl_v2 (more training samples, new genes, FOXO1 N-terminal feature; see
#' README_models.md of the GAMBL LymphGen repository).
#' @param not_assessed Optional gene symbols (or Entrez IDs) the samples were
#' not assessed for, e.g. not covered by their capture
#' (`Run_Lymphgen.R --not_assessed`, GAMBL LymphGen repository): left out of
#' the gene list so their features count as unknown rather than unmutated.
#' It applies to every sample of the run, so run cohorts with different
#' captures separately. For the Reddy et al. exomes:
#' `c("LTB", "HLA-A", "HLA-B", "P2RY8", "PRRC2A")`.
#' @param version LymphGen algorithm version: 2 (default; improved for
#' incomplete data) or 1 (the original).
#' @param a53 Classify cases into A53. Default TRUE.
#' @param cnv_mode Optional `Run_Lymphgen.R` copy number mode: "Full",
#' "None", "Only_AMP_HOMDEL" or "Only_Gain_Loss". Default NULL: "Full" with
#' copy number, "None" without.
#' @param compare Also return LymphGen's per-class precision and recall for
#' the features provided (`--compareout`), as `attr(result, "compare")`.
#' Default FALSE.
#' @param rscript Rscript to run LymphGen with. Default: the
#' `GAMBLR.predict.lymphgen_rscript` option, else this R's Rscript.
#'
#' @return For one input, a data frame with one row per sample: sample_id,
#' LymphGen (the subtype prediction, e.g. "EZB", "EZB/ST2" or "Other"), then
#' LymphGen's own columns (Copy.Number, the translocations, Model.Used, the
#' confidence, feature count and features of each subtype, and
#' Subtype.Prediction). For a list, a list of them with the same names. The directory, the
#' checksum of `Source LymphGen.R`, the options and the console output are in
#' `attr(result, "lymphgen_run")`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' options(GAMBLR.predict.lymphgen_dir = "~/git/Lymphgen_2")
#' inputs = lymphgen_input(maf, meta, seg_data = segs)
#' calls = run_lymphgen(inputs)
#' table(calls$LymphGen)
#'
#' # with and without copy number
#' no_cn = lymphgen_input(maf, meta)
#' calls = run_lymphgen(list(with_cn = inputs, no_cn = no_cn))
#'
#' # the GAMBL-derived model (see ?lymphgen_input for its inputs)
#' v2 = lymphgen_input(maf, meta, noncoding_maf = noncoding, model = "gambl_v2")
#' calls_v2 = run_lymphgen(v2)
#' plot_classification_alluvial(calls, class_column = "LymphGen")  # GAMBLR.viz
#' }
run_lymphgen = function(lymphgen_input,
                        lymphgen_dir = NULL,
                        model = NULL,
                        version = 2,
                        a53 = TRUE,
                        cnv_mode = NULL,
                        not_assessed = NULL,
                        compare = FALSE,
                        rscript = NULL){
  if(is.null(lymphgen_dir)){
    lymphgen_dir = getOption("GAMBLR.predict.lymphgen_dir", Sys.getenv("LYMPHGEN_DIR"))
  }
  if(is.null(rscript)){
    rscript = getOption("GAMBLR.predict.lymphgen_rscript", file.path(R.home("bin"), "Rscript"))
  }
  lymphgen_dir = path.expand(lymphgen_dir)
  rscript = path.expand(rscript)
  if(!version %in% c(1, 2)) stop("version must be 1 or 2")
  if(!is.null(cnv_mode)){
    cnv_mode = match.arg(cnv_mode, c("Full", "None", "Only_AMP_HOMDEL", "Only_Gain_Loss"))
  }
  needed = c("Run_Lymphgen.R", "Source LymphGen.R", "Fullmat.2020.txt", "Full.Uber.2020.txt",
             "Featord.2020.txt", "Index.Flat.txt", "Study.samp.txt", "Result.txt",
             if(version == 1) "originalpred.txt" else "originalpredv2.txt")
  if(identical(lymphgen_dir, "") || !dir.exists(lymphgen_dir)){
    stop("LymphGen directory not found", if(!identical(lymphgen_dir, "")) paste0(": ", lymphgen_dir) else "",
         "; give lymphgen_dir (or set options(GAMBLR.predict.lymphgen_dir = ...)); see ?run_lymphgen")
  }
  # the model: given, else the inputs', else the original
  single = is.list(lymphgen_input) &&
    any(names(lymphgen_input) %in% c("mutation_flat", "gene_list", "sample_annotation", "copy_flat"))
  input_models = unique(unlist(lapply(if(single) list(lymphgen_input) else lymphgen_input,
                                      function(x) attr(x, "lymphgen_model"))))
  if(is.null(model)){
    if(length(input_models) > 1) stop("the inputs were built for different models (", paste(input_models, collapse = ", "), "): give model")
    model = if(length(input_models) == 1) input_models else "original"
  }
  model = match.arg(model, c("original", "gambl_v2", "gambl_v3"))
  if(length(input_models) > 0 && !all(input_models == model)){
    warning("inputs built for model ", paste(input_models, collapse = ", "), " are run with model ", model, call. = FALSE)
  }
  if(model != "original"){
    needed = c(needed, file.path("models", model, "config.txt"))
    if(version != 2) stop("model ", model, " supports version 2 only")
  }
  absent = needed[!file.exists(file.path(lymphgen_dir, needed))]
  if(length(absent) > 0){
    stop(lymphgen_dir, " lacks LymphGen file(s): ", paste(absent, collapse = ", "),
         if(model != "original") " (models other than the original need the GAMBL LymphGen repository)" else "", "; see ?run_lymphgen")
  }
  has_optparse = suppressWarnings(system2(rscript, c("-e", shQuote("cat(requireNamespace('optparse', quietly = TRUE))")),
                                          stdout = TRUE, stderr = TRUE))
  if(!identical(utils::tail(has_optparse, 1), "TRUE")){
    stop("Run_Lymphgen.R needs the optparse package in ", rscript, ": install.packages(\"optparse\")")
  }

  # one input (lymphgen_input()'s list of data frames) or a named list of them
  inputs = if(single) list(lymphgen = lymphgen_input) else lymphgen_input
  # unnamed ones are named input1, input2, ...
  if(is.null(names(inputs))) names(inputs) = rep("", length(inputs))
  unnamed = names(inputs) == ""
  names(inputs)[unnamed] = paste0("input", which(unnamed))
  for(n in names(inputs)){
    parts = setdiff(c("mutation_flat", "gene_list", "sample_annotation"), names(inputs[[n]]))
    if(length(parts) > 0){
      stop(if(single) "lymphgen_input" else paste0("input ", n), " lacks ", paste(parts, collapse = ", "),
           " (use the output of lymphgen_input())")
    }
  }

  dir = tempfile("lymphgen_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  script = file.path(lymphgen_dir, "Run_Lymphgen.R")
  source_md5 = unname(tools::md5sum(file.path(lymphgen_dir, "Source LymphGen.R")))
  calls = lapply(seq_along(inputs), function(i){
    x = inputs[[i]]
    write_part = function(part){
      if(is.null(x[[part]])) return(NULL)
      path = file.path(dir, paste0(i, "_", part, ".tsv"))
      utils::write.table(x[[part]], path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
      path
    }
    files = lapply(c(mutation_flat = "mutation_flat", gene_list = "gene_list",
                     sample_annotation = "sample_annotation", copy_flat = "copy_flat", arm = "arm"),
                   write_part)
    out_file = file.path(dir, paste0(i, "_result.tsv"))
    compare_file = file.path(dir, paste0(i, "_compare.tsv"))
    args = c(shQuote(script),
             "-m", shQuote(files$mutation_flat), "-g", shQuote(files$gene_list),
             "-s", shQuote(files$sample_annotation),
             if(!is.null(files$copy_flat)) c("-c", shQuote(files$copy_flat)),
             if(!is.null(files$copy_flat) && !is.null(files$arm)) c("-a", shQuote(files$arm)),
             "-o", shQuote(out_file), "-v", version,
             if(model != "original") c("--model", model),
             if(!a53) "--no_A53",
             if(!is.null(cnv_mode)) c("--cnv_mode", cnv_mode),
             if(length(not_assessed) > 0) c("--not_assessed", shQuote(paste(not_assessed, collapse = ","))),
             if(compare) c("--compareout", shQuote(compare_file)))
    log = suppressWarnings(system2(rscript, args, stdout = TRUE, stderr = TRUE))
    if(!is.null(attr(log, "status")) || !file.exists(out_file)){
      stop("LymphGen failed", if(!single) paste0(" for ", names(inputs)[i]), ":\n",
           paste(utils::tail(log, 10), collapse = "\n"))
    }
    result = utils::read.delim(out_file, check.names = FALSE, stringsAsFactors = FALSE)
    result = data.frame(sample_id = result$Sample.Name, LymphGen = result$Subtype.Prediction,
                        result[, setdiff(colnames(result), "Sample.Name"), drop = FALSE],
                        check.names = FALSE, stringsAsFactors = FALSE)
    attr(result, "lymphgen_run") = list(lymphgen_dir = lymphgen_dir, source_md5 = source_md5,
                                        model = model, version = version, a53 = a53 && model == "original",
                                        cnv_mode = if(is.null(cnv_mode)) {
                                          if(is.null(files$copy_flat)) "None" else "Full"
                                        } else cnv_mode,
                                        log = log)
    if(compare && file.exists(compare_file)){
      attr(result, "compare") = utils::read.delim(compare_file, check.names = FALSE, stringsAsFactors = FALSE)
    }
    result
  })
  names(calls) = names(inputs)
  if(single) calls[[1]] else calls
}
