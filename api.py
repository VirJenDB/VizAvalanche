import httpx, logging, os.path, traceback

import pandas as pd
from pandas import DataFrame
from operator import itemgetter

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, JSONResponse

from config.config import config

import numpy as np
from scipy.stats import gaussian_kde
import plotly.graph_objs as go

import constants.fields as field
from constants.misc import HITS
from constants.misc   import GET
from shared_objects.elastic import es_async

LOG = logging.getLogger('GraphApi')
CUSTOM_COLORS = ["#E37B40", "#58A4B0"] * 20

REQUIRED_FIELDS = [
    field.CLUSTER_REPRESENTATIVE,
    field.MOLECULE_TYPE,
    field.NCBI_FAMILY,
    field.NCBI_GENUS,
    field.NCBI_SPECIES,
    field.SEQUENCE_GC_CONTENT,
    field.SEQUENCE_LENGTH,
    field.SEQUENCE_TITLE,
    field.SEQUENCE_COMPLETENESS,
    field.SUBMITTER_COUNTRY
]

def _base_layout(fig, xtitle, ytitle, tick_dist):
    fig.update_layout(
        #        height=500,
        #        width=700,
        bargap=.25,
        font=dict(size=16,color='black'),
        margin=dict(pad=16,t=40,b=10,l=0,r=0),
        plot_bgcolor='white',
        showlegend=False,
        xaxis=dict(
            title=dict(text='<b>{0}</b>'.format(xtitle)),
            showgrid=True,
            ticks="outside",
            tickson="boundaries",
            ticklen=5,
            tick0 = 0,
            nticks=7,
            #            dtick=tick_dist,
            showline=True,
            linewidth=1,
            linecolor='black',
            type='log',
            gridwidth=1,
            gridcolor='#cccccc'
        ),
        yaxis=dict(
            #            title=dict(
            #                text='<b>{0}</b>'.format(ytitle),
            #            ),
            showgrid=True,
            ticks="outside",
            tickson="labels",
            ticklen=5,
            showline=True,
            linewidth=1,
            linecolor='black',
        ),
    )

def _break_names(list_of_tuples):
    renamed_tuples = []
    for tuple in list_of_tuples:
        name = tuple[0]
        mid = len(name) // 2
        if mid < 7:
            # name already short
            tup = ('<sup>'+name+'</sup>',tuple[1])
            renamed_tuples.append(tup)
            continue
        left = name.rfind(' ', 0, mid+1)
        right = name.find(' ', mid)
        if left == -1 and right == -1:
            # No spaces found
            tup = ('<sub>'+name+'</sub>',tuple[1])
            renamed_tuples.append(tup)
            continue

        if left == -1 or (right != -1 and (right - mid) < (mid - left)):
            # Only right or closer right
            split_idx = right
        else:
            split_idx = left
        tup = ('<sub>'+name[:split_idx]+'</sub><br><sup>'+name[split_idx+1:]+'</sup>',tuple[1])
        renamed_tuples.append(tup)
    return renamed_tuples

def _create_density_gc_content(data_frame):
    gc = data_frame[field.SEQUENCE_GC_CONTENT][1:-1:10]
    select = np.logical_and(gc<100.0, gc>0.0)
    df_reps = data_frame[data_frame[field.CLUSTER_REPRESENTATIVE] == True]
    gc2 = df_reps[field.SEQUENCE_GC_CONTENT]
    select2 = np.logical_and(gc2<100.0, gc2>0.0)

    kde1 = gaussian_kde(select)
    kde2 = gaussian_kde(select2)
    x_grid1 = np.linspace(min(select), max(select), 1000)
    x_grid2 = np.linspace(min(select2), max(select2), 1000)
    density1 = kde1(x_grid1)
    density2 = kde2(x_grid2)

    fig = go.Figure()
    fig.add_trace(
        go.Line(x=x_grid1,y=density1)
    )
    fig.add_trace(
        go.Line(x=x_grid2,y=density2)
    )
    fig.update_layout(
        #        height=500,
        #        width=700,
        bargap=.25,
        font=dict(size=16,color='black'),
        margin=dict(pad=16,t=40,b=10,l=0,r=0),
        plot_bgcolor='white',
        showlegend=False,
        xaxis=dict(
            title=dict(text='<b>{0}</b>'.format('bla')),
            showgrid=True,
            ticks="outside",
            tickson="boundaries",
            ticklen=5,
            tick0 = 0,
            nticks=7,
            #            dtick=tick_dist,
            showline=True,
            linewidth=1,
            linecolor='black',
            type='log',
            gridwidth=1,
            gridcolor='#cccccc'
        ),
        yaxis=dict(
            #            title=dict(
            #                text='<b>{0}</b>'.format(ytitle),
            #            ),
            showgrid=True,
            ticks="outside",
            tickson="labels",
            ticklen=5,
            showline=True,
            linewidth=1,
            linecolor='black',
        ),
    )
    #fig.show()
    data_cache = f"{config.storage.plot_data}/density_gc_content.json"
    with open(data_cache, "w") as f:
        f.write(fig.to_json())
    LOG.debug(f'Plot {data_cache} created.')

def _create_density_seq_length(data_frame):
    lens = data_frame[field.SEQUENCE_LENGTH][1:-1:10]
    lens = lens[lens < 250000]
    df_reps = data_frame[data_frame[field.CLUSTER_REPRESENTATIVE] == True]
    lens2 = df_reps[field.SEQUENCE_LENGTH][df_reps[field.SEQUENCE_LENGTH] < 250000]

    kde1 = gaussian_kde(lens)
    kde2 = gaussian_kde(lens2)
    x_grid1 = np.linspace(min(lens), max(lens), 1000)
    x_grid2 = np.linspace(min(lens2), max(lens2), 1000)
    density1 = kde1(x_grid1)
    density2 = kde2(x_grid2)

    fig = go.Figure()
    fig.add_trace(
        go.Line(x=x_grid1,y=density1)
    )
    fig.add_trace(
        go.Line(x=x_grid2,y=density2)
    )
    fig.update_layout(
        #        height=500,
        #        width=700,
        bargap=.25,
        font=dict(size=16,color='black'),
        margin=dict(pad=16,t=40,b=10,l=0,r=0),
        plot_bgcolor='white',
        showlegend=False,
        xaxis=dict(
            title=dict(text='<b>{0}</b>'.format('bla')),
            showgrid=True,
            ticks="outside",
            tickson="boundaries",
            ticklen=5,
            tick0 = 0,
            nticks=7,
            #            dtick=tick_dist,
            showline=True,
            linewidth=1,
            linecolor='black',
            gridwidth=1,
            gridcolor='#cccccc'
        ),
        yaxis=dict(
            #            title=dict(
            #                text='<b>{0}</b>'.format(ytitle),
            #            ),
            showgrid=True,
            ticks="outside",
            tickson="labels",
            ticklen=5,
            showline=True,
            linewidth=1,
            linecolor='black',
            type='log',
        ),
    )
    #fig.show()

    data_cache = f"{config.storage.plot_data}/density_seq_length.json"
    with open(data_cache, "w") as f:
        f.write(fig.to_json())
    LOG.debug(f'Plot {data_cache} created.')

def _create_graph_country_stats(data_frame):
    data_cache = f"{config.storage.plot_data}/country_stats.json"
    max_bars       = 10
    threshold      = 0
    country_counts = data_frame[field.SUBMITTER_COUNTRY].value_counts()
    largest = country_counts[country_counts.values > threshold]
    while largest.size > max_bars:
        threshold = threshold +1
        largest = country_counts[country_counts.values > threshold]
    remaining_sum = country_counts[country_counts.values <= threshold].sum()


    country_counts_sorted       = sorted(largest.items(),key=itemgetter(1),reverse=False)
    country_counts_sorted.insert(0,('other',remaining_sum))

    keys, values = zip(*country_counts_sorted)
    bar_chart = go.Bar(x=values,y=keys,orientation='h',marker_color=CUSTOM_COLORS)

    fig = go.Figure([bar_chart])
    _base_layout(fig,'Count (log scale)','Top 10 Submission Countries',None)
    with open(data_cache, "w") as f:
        f.write(fig.to_json())
    LOG.debug(f'Plot {data_cache} created.')

def _create_graph_molecule_type(data_frame):
    data_cache = f"{config.storage.plot_data}/molecule_type.json"
    molecule_type_counts = data_frame[field.MOLECULE_TYPE].value_counts()

    unknown_sum = molecule_type_counts[molecule_type_counts.index == 'unknown'].sum()
    molecule_type_counts = molecule_type_counts[molecule_type_counts.index != 'unknown']
    plot_data = sorted(molecule_type_counts.items(),key=itemgetter(1),reverse=False)
    plot_data.insert(0,('unknown',unknown_sum))
    keys, values = zip(*plot_data)
    bar_chart = go.Bar(x=values, y=keys, orientation='h', marker_color=CUSTOM_COLORS)

    fig = go.Figure([bar_chart])
    _base_layout(fig, 'Count (log scale)', 'Molecule Type', None)
    with open(data_cache, "w") as f:
        f.write(fig.to_json())
    LOG.debug(f'Plot {data_cache} created.')

def _create_graph_sample_family(data_frame):
    data_cache = f"{config.storage.plot_data}/sample_family.json"
    max_bars = 10
    threshold = 0
    family_counts_inter  = data_frame[field.NCBI_FAMILY].value_counts()

    largest = family_counts_inter[family_counts_inter.values > threshold]
    while largest.size > max_bars:
        threshold += 1
        largest = family_counts_inter[family_counts_inter.values > threshold]
    remaining_sum = family_counts_inter[family_counts_inter.values <= threshold].sum()
    family_counts        = sorted(
        largest.items(),
        key=lambda x: (-x[1], x[0]), # sort by frequency, then name
        reverse=True) # reverse order
    family_counts.insert(0,('other',remaining_sum))
    keys, values = zip(*family_counts)
    bar_chart = go.Bar(x=values,y=keys,orientation='h',marker_color=CUSTOM_COLORS)

    fig = go.Figure([bar_chart])
    _base_layout(fig,'Count (log scale)','Top 10 NCBI Families',None)
    fig.update_traces(marker_color=CUSTOM_COLORS)

    with open(data_cache, "w") as f:
        f.write(fig.to_json())
    LOG.debug(f'Plot {data_cache} created.')

def _create_graph_sample_genus(data_frame):
    data_cache = f"{config.storage.plot_data}/sample_genus.json"
    max_bars = 10
    threshold = 0
    genus_counts_inter  = data_frame[field.NCBI_GENUS].value_counts()

    largest = genus_counts_inter[genus_counts_inter.values > threshold]
    while largest.size > max_bars:
        threshold += 1
        largest = genus_counts_inter[genus_counts_inter.values > threshold]
    remaining_sum = genus_counts_inter[genus_counts_inter.values <= threshold].sum()
    genus_counts        = sorted(
        largest.items(),
        key=lambda x: (-x[1], x[0]), # sort by frequency, then name
        reverse=True) # reverse order
    genus_counts.insert(0,('other',remaining_sum))
    keys, values = zip(*genus_counts)
    bar_chart = go.Bar(x=values,y=keys,orientation='h',marker_color=CUSTOM_COLORS)

    fig = go.Figure([bar_chart])
    _base_layout(fig,'Count (log scale)','Top 10 NCBI Genus',None)
    fig.update_traces(marker_color=CUSTOM_COLORS)

    with open(data_cache, "w") as f:
        f.write(fig.to_json())
    LOG.debug(f'Plot {data_cache} created.')

def _create_graph_sample_species(data_frame):
    data_cache = f"{config.storage.plot_data}/sample_species.json"
    max_bars = 10
    threshold = 0
    species_counts_inter  = data_frame[field.NCBI_SPECIES].value_counts()

    largest = species_counts_inter[species_counts_inter.values > threshold]
    while largest.size > max_bars:
        threshold += 1
        largest = species_counts_inter[species_counts_inter.values > threshold]
    remaining_sum = species_counts_inter[species_counts_inter.values <= threshold].sum()
    species_counts        = sorted(
        largest.items(),
        key=lambda x: (-x[1], x[0]), # sort by frequency, then name
        reverse=True) # reverse order
    species_counts.insert(0,('other',remaining_sum))

    renamed_species_counts = _break_names(species_counts)

    keys, values = zip(*renamed_species_counts)
    bar_chart = go.Bar(x=values,y=keys,orientation='h',marker_color=CUSTOM_COLORS)

    fig = go.Figure([bar_chart])
    _base_layout(fig,'Count (log scale)','Top 10 NCBI species',None)
    fig.update_traces(marker_color=CUSTOM_COLORS)

    with open(data_cache, "w") as f:
        f.write(fig.to_json())
    LOG.debug(f'Plot {data_cache} created.')

def _create_graph_seq_len_by_family(data_frame):
    data_cache = f"{config.storage.plot_data}/seq_len_by_family.json"
    data_frame = data_frame[data_frame[field.SEQUENCE_COMPLETENESS] == 'complete']
    order = (
        data_frame.groupby(field.NCBI_FAMILY)[field.SEQUENCE_LENGTH]
        .median()
        .sort_values()
        .index
    )
    fig = go.Figure()
    for family in order:
        family_data = data_frame[data_frame[field.NCBI_FAMILY] == family]
        seq_len = family_data[field.SEQUENCE_LENGTH].to_list()
        fig.add_trace(
            go.Box(
                y=seq_len,
                name=family,
                boxpoints='outliers', # only outliers
                marker=dict(color='black'),
                fillcolor='white',
                line=dict(color='black'),
                jitter=0,
            )
        )

    _base_layout(fig,'Family','Sequence Length',None)
    fig.update_layout(
        xaxis=dict(tickangle=-45,tickson='labels',type='linear'),
        yaxis=dict(
            type="log",
            #exponentformat="e",     # Show exponent as "e" (e.g. 1e-5 not 10^-5)
            #showexponent="all",      # Always show exponent, even if not strictly required
            #range=[-2e5,4e5]
        ),
    )
    with open(data_cache, "w") as f:
        f.write(fig.to_json())

async def _fill_data_cache(csv_file) -> DataFrame:
    cache_time   = "10m"
    requests     = 0
    no_seq_title = 0
    try:
        query = { "match_all" : {} }
        collected = []
        response  = await es_async.search(index=config.elastic.index.main,scroll=cache_time,source=REQUIRED_FIELDS,size=10000,query=query)
        while True:
            data      = response.body
            requests += 1
            if not (isinstance(data.get(HITS),dict) and isinstance(data[HITS].get(HITS),list) and data[HITS][HITS]):
                break

            hits  = data[HITS][HITS]

            for hit in hits: # go through entries of returned data slice
                interest  = []
                for idx,field_id in enumerate(REQUIRED_FIELDS): # construct dict of relevant fields for each entry
                    value = hit['_source'][field_id] if hit['_source'].get(field_id) else None
                    if value is None:
                        if field_id == field.SEQUENCE_TITLE:
                            no_seq_title +=1
                            break
                    else: # value is present
                        if field_id == field.SEQUENCE_LENGTH:
                            value = int(value)
                        if field_id == field.SEQUENCE_GC_CONTENT:
                            value = float(value)
                    interest.append(value)
                else: # interest contains no sequence title, so don't add it to collection
                    collected.append(interest) # collect data

            LOG.debug(f'fired {requests} requests, received {len(collected)} entries with sequence title ({no_seq_title} entries have no sequence_title)')

            response = await es_async.scroll(scroll_id=data['_scroll_id'],scroll=cache_time)
        data_frame = pd.DataFrame(columns=REQUIRED_FIELDS,data=collected)
        data_frame.to_csv(csv_file)
        return data_frame
    except httpx.HTTPStatusError as exc:
        raise HTTPException(status_code=exc.response.status_code, detail=exc.response.json())
    except HTTPException as hx:
        raise hx
    except Exception as e:
        print(traceback.format_exc())
        raise HTTPException(status_code=500, detail=f"Unexpected error: {str(e)}")

async def _get_graph_molecule_type():
    data_cache = f"{config.storage.plot_data}/molecule_type.json"
    if os.path.isfile(data_cache):
        return FileResponse(path=data_cache,media_type='application/json')
    raise HTTPException(status_code=503, detail='Plot data needs recreation!')



async def _get_graph_sample_family():
    data_cache = f"{config.storage.plot_data}/sample_family.json"
    if os.path.isfile(data_cache):
        return FileResponse(path=data_cache,media_type='application/json')
    raise HTTPException(status_code=503, detail='Plot data needs recreation!')



async def _get_graph_sample_genus():
    data_cache = f"{config.storage.plot_data}/sample_genus.json"
    if os.path.isfile(data_cache):
        return FileResponse(path=data_cache,media_type='application/json')
    raise HTTPException(status_code=503, detail='Plot data needs recreation!')



async def _get_graph_sample_species():
    data_cache = f"{config.storage.plot_data}/sample_species.json"
    if os.path.isfile(data_cache):
        return FileResponse(path=data_cache,media_type='application/json')
    raise HTTPException(status_code=503, detail='Plot data needs recreation!')



async def _get_graph_country_stats():
    data_cache = f"{config.storage.plot_data}/country_stats.json"
    if os.path.isfile(data_cache):
        return FileResponse(path=data_cache,media_type='application/json')
    raise HTTPException(status_code=503, detail='Plot data needs recreation!')



async def _get_graph_seq_len_by_family():
    data_cache = f"{config.storage.plot_data}/seq_len_by_family.json"
    if os.path.isfile(data_cache):
        return FileResponse(path=data_cache,media_type='application/json')
    raise HTTPException(status_code=503, detail='Plot data needs recreation!')


async def _get_graph_density_seq_length():
    data_cache = f"{config.storage.plot_data}/density_seq_length.json"
    if os.path.isfile(data_cache):
        return FileResponse(path=data_cache,media_type='application/json')
    raise HTTPException(status_code=503, detail='Plot data needs recreation!')


async def _get_graph_density_gc_content():
    data_cache = f"{config.storage.plot_data}/density_gc_content.json"
    if os.path.isfile(data_cache):
        return FileResponse(path=data_cache,media_type='application/json')
    raise HTTPException(status_code=503, detail='Plot data needs recreation!')


async def _get_plot_data():
    data_cache = f"{config.storage.plot_data}/collected_data.csv"
    if os.path.isfile(data_cache):
        LOG.debug(f'reading data from {data_cache}. This may take a while…')
        df = pd.read_csv(data_cache)
        LOG.debug(f'…done')
        return df
    return await _fill_data_cache(data_cache)


async def _recreate_plot_data():
    data_frame = await _get_plot_data()
    _create_graph_molecule_type(data_frame)
    _create_graph_sample_family(data_frame)
    _create_graph_sample_genus(data_frame)
    _create_graph_country_stats(data_frame)
    _create_graph_sample_species(data_frame)
    _create_density_seq_length(data_frame)
    _create_density_gc_content(data_frame)
    #_create_graph_seq_len_by_family(data_frame)
    return "done"


class GraphApi(APIRouter):
    def __init__(self):
        super().__init__()
        self.add_api_route(path="/recreate",      endpoint=_recreate_plot_data,      methods=[GET], response_class=JSONResponse, tags=['Graphs'])

        common_args = {
            'methods' : [GET],
            'response_class' : JSONResponse,
            'tags' : ['Graphs']
        }
        self.add_api_route(path="/molecule_type",      endpoint=_get_graph_molecule_type,      **common_args)
        self.add_api_route(path="/sample_family",      endpoint=_get_graph_sample_family,      **common_args)
        self.add_api_route(path="/sample_genus",       endpoint=_get_graph_sample_genus,       **common_args)
        self.add_api_route(path="/sample_species",     endpoint=_get_graph_sample_species,     **common_args)
        self.add_api_route(path="/country_stats",      endpoint=_get_graph_country_stats,      **common_args)
        self.add_api_route(path="/seq_len_by_family",  endpoint=_get_graph_seq_len_by_family,  **common_args)
        self.add_api_route(path="/density_seq_length", endpoint=_get_graph_density_seq_length, **common_args)
        self.add_api_route(path="/density_gc_content", endpoint=_get_graph_density_gc_content, **common_args)







