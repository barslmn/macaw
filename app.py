#!/usr/bin/env python3
from collections import defaultdict
from io import BytesIO, StringIO
from Bio import SeqIO
from Bio import Align

import base64
from dash import Dash, html, dcc, Input, Output
import plotly.graph_objs as go
from plotly.subplots import make_subplots

index_offset = 200


def get_theme(theme):
    themes = {
        "default": {
            "black": "#000000",
            "red": "#ff0000",
            "green": "#00ff00",
            "blue": "#0000ff",
            "orange": "#ff8000",
            "cyan": "#00ffff",
            "violet": "#8000ff",
            "pink": "#ff0080",
        },
        "excalidraw": {
            "black": "#1e1e1e",
            "red": "#e03131",
            "green": "#2f9e44",
            "blue": "#1971c2",
            "orange": "#f08c00",
            "cyan": "#0c8599",
            "violet": "#6741d9",
            "pink": "#c2255c",
        },
    }
    return themes[theme]


chromatogram_theme = get_theme("excalidraw")


def align(seq1, seq2):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    return alignments[0]


def find_cealing(channels, trace):
    cealing = 0
    for channel in channels:
        if max(trace[channel]) > cealing:
            cealing = max(trace[channel])
    return cealing


def add_indel_peaks(peaks, trace_seq):
    for i, char in enumerate(trace_seq):
        # slide every peak to right for every deletion
        if char == "-":
            peak = peaks[i]
            if i > 0:
                peak_diff = peak - peaks[i - 1]
            else:
                peak_diff = peak
            peaks = peaks[:i] + [p + peak_diff for p in peaks[i:]]
            peaks.insert(i, peak)
    return peaks


def mark_variation(variation, peaks, cealing, index):
    if not variation:
        return []
    var_marks = []
    for i, char in enumerate(variation):
        x0 = peaks[i] - 10
        y0 = 0
        if char == " ":
            fillcolor = "rgba(255, 0, 0, 0.2)"
            color = "rgba(255, 0, 0, 0.2)"
        elif char == ".":
            fillcolor = "rgba(255, 165, 0, 0.2)"
            color = "rgba(255, 165, 0, 0.2)"
        else:
            continue

        # TO DO better fit peaks
        var_marks.append(
            {
                "type": "rect",
                "xref": "x{}".format(index + 1),
                "yref": "y{}".format(index + 1),
                "x0": x0,
                "y0": y0,
                "x1": x0 + 20,
                "y1": y0 + cealing + index_offset,
                "line": {
                    "color": color,
                    "width": 2,
                },
                "fillcolor": fillcolor,
            }
        )
    return var_marks


def widen_trace(data):
    return [i for i in range(len(data))]


def process_trace(data, trace_seq, peaks):
    data = list(data)
    for i, char in enumerate(trace_seq):
        if char == "-":
            data[peaks[i] : peaks[i]] = [0] * (peaks[i + 1] - peaks[i])
    return data[: peaks[-1]]


def preprocess_chromatogram(abi_seq, ref_seq=None):
    record = SeqIO.read(BytesIO(abi_seq), "abi")

    peaks = list(record.annotations["abif_raw"]["PLOC1"])

    channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
    trace = defaultdict(list)

    for channel in channels:
        trace[channel] = record.annotations["abif_raw"][channel]

    cealing = find_cealing(channels, trace)

    variation, trace_seq = None, None
    if ref_seq is not None:
        ref_seq = SeqIO.read(StringIO(ref_seq), "fasta").seq.upper()
        alignments = align(ref_seq, record.seq)
        variation = "".join(
            map(
                lambda x: x[1],
                filter(
                    lambda x: len(x) > 1,
                    map(
                        lambda x: x.strip().split(),
                        filter(None, alignments[0].format().split("\n")[1::2]),
                    ),
                ),
            )
        )
        ref_seq, trace_seq = alignments

        peaks = add_indel_peaks(peaks, trace_seq)

    return {
        "record": record,
        "trace": trace,
        "trace_seq": trace_seq,
        "variation": variation,
        "ref_seq": ref_seq,
        "peaks": peaks,
        "cealing": cealing,
    }


def get_nucleotide_color(nucleotide):
    colors = {
        "A": chromatogram_theme["green"],
        "T": chromatogram_theme["red"],
        "C": chromatogram_theme["blue"],
        "G": chromatogram_theme["black"],
        "-": chromatogram_theme["black"],
    }
    return f"<span style='color:{colors[nucleotide]}'> {nucleotide} </span>"


def create_chromatogram(abi_contents, ref_content={"ref_content": None}):

    subplots = []
    for i, (abi_name, abi_content) in enumerate(abi_contents.items()):

        chromatogram_data = preprocess_chromatogram(
            abi_content, ref_content["ref_content"]
        )
        trace = chromatogram_data.get("trace", None)
        record = chromatogram_data.get("record", None)
        peaks = chromatogram_data.get("peaks", None)
        cealing = chromatogram_data.get("cealing", None)

        trace_seq = chromatogram_data.get("trace_seq", None)
        variation = chromatogram_data.get("variation", None)
        ref_seq = chromatogram_data.get("ref_seq", None)

        if not trace_seq:
            trace_seq = record.seq

        trace1 = go.Scatter(
            x=widen_trace(trace["DATA9"]),
            y=process_trace(trace["DATA9"], trace_seq, peaks),
            name="{}".format(chr(record.annotations["abif_raw"]["FWO_1"][0])),
            line_shape="spline",
            marker=dict(color=chromatogram_theme["black"]),
        )
        trace2 = go.Scatter(
            x=widen_trace(trace["DATA10"]),
            y=process_trace(trace["DATA10"], trace_seq, peaks),
            name="{}".format(chr(record.annotations["abif_raw"]["FWO_1"][1])),
            line_shape="spline",
            marker=dict(color=chromatogram_theme["green"]),
        )
        trace3 = go.Scatter(
            x=widen_trace(trace["DATA11"]),
            y=process_trace(trace["DATA11"], trace_seq, peaks),
            name="{}".format(chr(record.annotations["abif_raw"]["FWO_1"][2])),
            line_shape="spline",
            marker=dict(color=chromatogram_theme["red"]),
        )
        trace4 = go.Scatter(
            x=widen_trace(trace["DATA12"]),
            y=process_trace(trace["DATA12"], trace_seq, peaks),
            name="{}".format(chr(record.annotations["abif_raw"]["FWO_1"][3])),
            line_shape="spline",
            marker=dict(color=chromatogram_theme["blue"]),
        )
        trace5 = go.Scatter(
            x=peaks[::10],
            y=[cealing + index_offset for i in range(len(peaks[: len(trace_seq)]))],
            mode="text",
            name="Index",
            text=[str(i) for i in range(1, len(peaks[: len(trace_seq)]), 10)],
            textposition="bottom center",
            marker=dict(color=plot_bgcolor),
            textfont=dict(color=chromatogram_theme["black"]),
        )

        data = [trace1, trace2, trace3, trace4, trace5]
        if ref_seq is not None:

            trace6 = go.Scatter(
                x=peaks[: len(ref_seq)],
                y=[cealing + 300 for i in range(len(peaks[: len(ref_seq)]))],
                # mode='markers+text',
                mode="text",
                name="Reference Sequence",
                text=list(ref_seq),
                textposition="bottom center",
                marker=dict(color=plot_bgcolor),
                textfont=dict(color=chromatogram_theme["orange"]),
            )
            data += [trace6]

        # info = [
        #     html.H2(children='{}'.format(record.name)),
        #     html.Div(children='''
        #         Machine Model: {} Run Start: {} Run Finish: {}
        #     '''.format(record.annotations['machine_model'], record.annotations['run_start'], record.annotations['run_finish'])),
        # ]
        subplot = dict(
            data=data,
            layout=go.Layout(
                plot_bgcolor="#ffffff",
                paper_bgcolor="#ffffff",
            ),
            xaxis=dict(
                # autorange=True,
                showgrid=False,
                rangeslider=dict(),
                # zeroline=False,
                # showline=False,
                # showticklabels=False
                tickmode="array",
                tickvals=peaks[: len(record.seq)],
                ticktext=[get_nucleotide_color(n) for n in trace_seq],
                tickangle=0,
                color="rgba(0, 0, 0, 1)",
            ),
            yaxis=dict(
                autorange=True,
                showgrid=False,
                zeroline=False,
                showline=False,
                tickmode="auto",
                ticks="",
                showticklabels=False,
                title=go.layout.yaxis.Title(
                    text=abi_name.split(".")[0],
                    font=dict(
                        # family='Courier New, monospace',
                        # size=18,
                        # color='#7f7f7f'
                    ),
                ),
            ),
            shapes=mark_variation(variation, peaks, cealing, i),
        )
        subplots.append(subplot)
    figure = make_subplots(rows=len(subplots), cols=1, shared_xaxes=True)
    # figure.update_yaxes(fixedrange=True)
    # figure.update_layout(yaxis_range=[0, cealing])
    for i, subplot in enumerate(subplots):
        for trace in subplot["data"]:
            figure.append_trace(trace, row=i + 1, col=1)
        figure.update_layout(subplot["layout"], shapes=subplot["shapes"])
        figure.update_xaxes(subplot["xaxis"])
        figure.update_yaxes(subplot["yaxis"])
    return figure


app = Dash(__name__)
plot_bgcolor = "#ffffff"


app.layout = html.Div(
    children=[
        html.Div(
            [
                html.Div(
                    [
                        html.H2("Abi file(s)"),
                        dcc.Upload(
                            id="upload-abi",
                            children=html.Div(
                                ["Drag and drop or click to select a file to upload."]
                            ),
                            style={
                                "width": "100%",
                                "height": "60px",
                                "lineHeight": "60px",
                                "borderWidth": "1px",
                                "borderStyle": "dashed",
                                "borderRadius": "5px",
                                "textAlign": "center",
                                "margin": "10px",
                            },
                            multiple=True,
                        ),
                    ],
                    className="six columns",
                ),
                html.Div(
                    [
                        html.H2("Reference file"),
                        dcc.Upload(
                            id="upload-ref",
                            children=html.Div(
                                ["Drag and drop or click to select a file to upload."]
                            ),
                            style={
                                "width": "100%",
                                "height": "60px",
                                "lineHeight": "60px",
                                "borderWidth": "1px",
                                "borderStyle": "dashed",
                                "borderRadius": "5px",
                                "textAlign": "center",
                                "margin": "10px",
                            },
                            multiple=False,
                        ),
                    ],
                    className="six columns",
                ),
                html.Div([dcc.Dropdown(["default", "excalidraw"], "default")]),
            ],
            className="row",
        ),
        # html.Div(id='chromatogram-info'),
        # html.Div(id='info'),
        html.Div(id="out"),
        # html.Div([
        #     dcc.Graph(id='graph'),
        # ])
    ]
)
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})


@app.callback(
    # Output("graph", "figure"),
    Output("out", "children"),
    [
        Input("upload-abi", "filename"),
        Input("upload-abi", "contents"),
        Input("upload-ref", "filename"),
        Input("upload-ref", "contents"),
    ],
)
def add_files(abi_names, abi_files, ref_name, ref_file):
    if abi_files or ref_file is not None:
        if ref_file:
            ref_content = {
                "ref_name": ref_name,
                "ref_content": base64.b64decode(ref_file.split(",")[1]).decode("ascii"),
            }
        else:
            ref_content = {"ref_content": None}
        if abi_files:
            abi_contents = {
                abi_name: base64.decodebytes(
                    abi_content.encode("utf8").split(b";base64,")[1]
                )
                for abi_name, abi_content in zip(abi_names, abi_files)
            }

        figure = create_chromatogram(abi_contents, ref_content=ref_content)
        return dcc.Graph(figure=figure)
    else:
        return [html.Div("No files yet!")]


if __name__ == "__main__":
    app.run(debug=True)
